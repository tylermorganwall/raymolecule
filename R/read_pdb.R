#' Read PDB File
#'
#' Reads a legacy fixed-width PDB file and extracts atom locations, bond
#' connections, backbone residue information, and secondary structure records.
#' `model$atoms` and `model$bonds` remain compatible with the existing atom and
#' bond scene generators.
#'
#' @param filename Path to the PDB file.
#' @param atom Default `TRUE`. Whether to include standard residue `ATOM`
#'   records in `model$atoms`.
#' @param nsr Default `TRUE`. Whether to include `HETATM` records in
#'   `model$atoms`.
#'
#' @return List giving the parsed PDB model.
#' @export
#'
#' @examples
#' # This assumes a hypothetical PDB file in your working directory:
#' if (file.exists("3nir.pdb")) {
#'   read_pdb("3nir.pdb") |>
#'     generate_full_scene() |>
#'     render_model()
#' }
read_pdb = function(filename, atom = TRUE, nsr = TRUE) {
	#Read in all the lines first
	lines = readLines(filename, warn = FALSE)
	#Parse the file
	records = parse_pdb_records(lines)

	atom_mask =
		(records$atoms$record == "ATOM" & atom) |
		(records$atoms$record == "HETATM" & nsr)

	atoms = records$atoms[
		atom_mask,
		c(
			"x",
			"y",
			"z",
			"type",
			"index",
			"record",
			"atom_name",
			"res_name",
			"chain_id",
			"res_seq",
			"i_code",
			"alt_loc"
		),
		drop = FALSE
	]

	return(list(
		atoms = atoms,
		bonds = records$bonds,
		residues = records$residues,
		helices = records$helices,
		sheets = records$sheets,
		ssbonds = records$ssbonds,
		chains = records$chains,
		pdb_type = "pdb"
	))
}

#' @keywords internal
parse_pdb_records = function(lines) {
	#Set up lists for all different objects within the ribbon to render
	atom_rows = list()
	bond_rows = list()
	helix_rows = list()
	sheet_rows = list()
	ssbond_rows = list()
	ter_rows = list()

	#Start counting
	atom_counter = 1L
	bond_counter = 1L
	helix_counter = 1L
	sheet_counter = 1L
	ssbond_counter = 1L
	ter_counter = 1L
	atom_order = 0L
	model_count = 0L
	in_model = FALSE

	# Parse each line
	for (line_number in seq_along(lines)) {
		line = lines[[line_number]]
		#Load the first 6 characters to get the record type
		record = trimws(substr(line, 1L, 6L))

		# Only support a single model (for now)
		if (identical(record, "MODEL")) {
			model_count = model_count + 1L
			if (model_count > 1L) {
				stop("Multiple MODEL records are not supported yet")
			}
			in_model = TRUE
			next
		}

		# We're done
		if (identical(record, "ENDMDL")) {
			in_model = FALSE
			next
		}

		# Are we in a model yet?
		active_model = model_count == 0L || in_model

		if (record %in% c("ATOM", "HETATM")) {
			# Don't parse if we aren't in a model
			if (!active_model) {
				next
			}
			atom_order = atom_order + 1L
			# Parse the atom record
			atom_rows[[atom_counter]] = parse_pdb_atom_record(
				line = line,
				line_number = line_number,
				record = record,
				atom_order = atom_order
			)
			atom_counter = atom_counter + 1L
			next
		}

		#Parsing bonds -- why don't we skip these below if we aren't active in a model as well?
		if (identical(record, "CONECT")) {
			bond_data = parse_pdb_conect_record(
				line = line,
				line_number = line_number
			)
			if (nrow(bond_data) > 0L) {
				bond_rows[[bond_counter]] = bond_data
				bond_counter = bond_counter + 1L
			}
			next
		}

		# PParse helix record
		if (identical(record, "HELIX")) {
			helix_rows[[helix_counter]] = parse_pdb_helix_record(
				line = line,
				line_number = line_number
			)
			helix_counter = helix_counter + 1L
			next
		}

		# Parse sheet record
		if (identical(record, "SHEET")) {
			sheet_rows[[sheet_counter]] = parse_pdb_sheet_record(
				line = line,
				line_number = line_number
			)
			sheet_counter = sheet_counter + 1L
			next
		}

		# Parse SSbond (?) record
		if (identical(record, "SSBOND")) {
			ssbond_rows[[ssbond_counter]] = parse_pdb_ssbond_record(
				line = line,
				line_number = line_number
			)
			ssbond_counter = ssbond_counter + 1L
			next
		}

		# Parse terminator record
		if (identical(record, "TER")) {
			if (!active_model) {
				next
			}
			ter_rows[[ter_counter]] = parse_pdb_ter_record(
				line = line,
				line_number = line_number
			)
			ter_counter = ter_counter + 1L
			next
		}
	}

	# This specifies the schema for each data type, and ensures it's enforced even at zero rows
	atoms = bind_pdb_rows(atom_rows, empty_pdb_atoms())
	bonds = bind_pdb_rows(bond_rows, empty_pdb_bonds())
	helices = bind_pdb_rows(helix_rows, empty_pdb_helices())
	sheets = bind_pdb_rows(sheet_rows, empty_pdb_sheets())
	ssbonds = bind_pdb_rows(ssbond_rows, empty_pdb_ssbonds())
	ters = bind_pdb_rows(ter_rows, empty_pdb_ters())

	# This builds the actual ribbon diagram data, with the helix, sheet, and terminator data
	residue_info = build_pdb_residues(
		atoms = atoms[atoms$record == "ATOM", , drop = FALSE],
		helices = helices,
		sheets = sheets,
		ters = ters
	)

	return(list(
		atoms = atoms,
		bonds = bonds,
		residues = residue_info$residues,
		chains = residue_info$chains,
		helices = helices,
		sheets = sheets,
		ssbonds = ssbonds
	))
}

#' @keywords internal
build_pdb_residues = function(atoms, helices, sheets, ters) {
	#Enforce schema
	residue_template = empty_pdb_residues()
	chain_template = empty_pdb_chains()

	if (nrow(atoms) == 0L) {
		return(list(residues = residue_template, chains = chain_template))
	}

	# Go though atoms
	residue_keys = paste(
		atoms$chain_id,
		atoms$res_seq,
		atoms$i_code,
		sep = "\r"
	)
	residue_key_order = unique(residue_keys)
	residue_rows = vector(mode = "list", length = length(residue_key_order))

	# Separate by residue, and extract N/CA/C/O atoms to get direction,
	# assign ss_class = "loop"
	for (i in seq_along(residue_key_order)) {
		residue_atoms = atoms[
			residue_keys == residue_key_order[[i]],
			,
			drop = FALSE
		]
		residue_atoms = residue_atoms[
			order(residue_atoms$atom_order),
			,
			drop = FALSE
		]
		first_atom = residue_atoms[1, , drop = FALSE]

		n_atom = select_residue_atom(residue_atoms, "N")
		ca_atom = select_residue_atom(residue_atoms, "CA")
		c_atom = select_residue_atom(residue_atoms, "C")
		o_atom = select_residue_atom(residue_atoms, "O")

		residue_rows[[i]] = data.frame(
			chain_id = first_atom$chain_id,
			res_seq = first_atom$res_seq,
			i_code = first_atom$i_code,
			res_name = first_atom$res_name,
			residue_id = format_pdb_residue_id(
				chain_id = first_atom$chain_id,
				res_seq = first_atom$res_seq,
				i_code = first_atom$i_code
			),
			n_x = if (is.null(n_atom)) NA_real_ else n_atom$x,
			n_y = if (is.null(n_atom)) NA_real_ else n_atom$y,
			n_z = if (is.null(n_atom)) NA_real_ else n_atom$z,
			ca_x = if (is.null(ca_atom)) NA_real_ else ca_atom$x,
			ca_y = if (is.null(ca_atom)) NA_real_ else ca_atom$y,
			ca_z = if (is.null(ca_atom)) NA_real_ else ca_atom$z,
			c_x = if (is.null(c_atom)) NA_real_ else c_atom$x,
			c_y = if (is.null(c_atom)) NA_real_ else c_atom$y,
			c_z = if (is.null(c_atom)) NA_real_ else c_atom$z,
			o_x = if (is.null(o_atom)) NA_real_ else o_atom$x,
			o_y = if (is.null(o_atom)) NA_real_ else o_atom$y,
			o_z = if (is.null(o_atom)) NA_real_ else o_atom$z,
			has_n = !is.null(n_atom),
			has_ca = !is.null(ca_atom),
			has_c = !is.null(c_atom),
			has_o = !is.null(o_atom),
			chain_break_before = FALSE,
			chain_break_after = FALSE,
			ss_class = "loop",
			helix_id = NA_character_,
			sheet_id = NA_character_,
			stringsAsFactors = FALSE
		)
	}

	residues = do.call(rbind, residue_rows)
	rownames(residues) = NULL

	residue_break_after = rep(FALSE, nrow(residues))
	residue_break_before = rep(FALSE, nrow(residues))
	ter_after = rep(FALSE, nrow(residues))

	if (nrow(residues) > 0L) {
		residue_break_before[1] = TRUE
		residue_break_after[nrow(residues)] = TRUE
	}

	if (nrow(residues) > 1L) {
		# Get when the chain changes
		chain_change = residues$chain_id[-1] !=
			residues$chain_id[-nrow(residues)]
		# Mark the ends
		residue_break_after[-nrow(residues)] =
			residue_break_after[-nrow(residues)] | chain_change
		# Mark the beginnings
		residue_break_before[-1] =
			residue_break_before[-1] | chain_change
	}

	if (nrow(ters) > 0L) {
		residue_match = match(
			paste(ters$chain_id, ters$res_seq, ters$i_code, sep = "\r"),
			paste(
				residues$chain_id,
				residues$res_seq,
				residues$i_code,
				sep = "\r"
			)
		)
		residue_match = residue_match[!is.na(residue_match)]
		if (length(residue_match) > 0L) {
			residue_break_after[residue_match] = TRUE
			ter_after[residue_match] = TRUE
			next_rows = residue_match[residue_match < nrow(residues)] + 1L
			if (length(next_rows) > 0L) {
				residue_break_before[next_rows] = TRUE
			}
		}
	}

	residues$chain_break_before = residue_break_before
	residues$chain_break_after = residue_break_after

	# Assign HELIX records
	residues = apply_pdb_secondary_structure(
		residues = residues,
		helices = helices
	)
	residues = apply_pdb_sheet_structure(residues = residues, sheets = sheets)

	chain_rows = list()
	chain_counter = 1L
	chain_start = 1L
	for (i in seq_len(nrow(residues))) {
		if (!residues$chain_break_after[i]) {
			next
		}

		chain_rows[[chain_counter]] = data.frame(
			chain_id = residues$chain_id[chain_start],
			chain_index = chain_counter,
			start_residue_id = residues$residue_id[chain_start],
			end_residue_id = residues$residue_id[i],
			terminated_by_ter = ter_after[i],
			stringsAsFactors = FALSE
		)
		chain_counter = chain_counter + 1L
		chain_start = i + 1L
	}
	# Extract chain start/end residue IDs
	chains = bind_pdb_rows(chain_rows, chain_template)

	return(list(residues = residues, chains = chains))
}

#' @keywords internal
parse_pdb_atom_record = function(line, line_number, record, atom_order) {
	atom_name = trimws(pdb_field(line, 13L, 16L, record, line_number))
	if (!nzchar(atom_name)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: blank atom name in columns 13-16",
				record,
				line_number
			)
		)
	}

	res_name = trimws(pdb_field(line, 18L, 20L, record, line_number))
	if (!nzchar(res_name)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: blank residue name in columns 18-20",
				record,
				line_number
			)
		)
	}

	element = trimws(pdb_field(line, 77L, 78L, record, line_number))
	if (!nzchar(element)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: blank element field in columns 77-78",
				record,
				line_number
			)
		)
	}

	return(data.frame(
		x = pdb_numeric_field(line, 31L, 38L, record, line_number),
		y = pdb_numeric_field(line, 39L, 46L, record, line_number),
		z = pdb_numeric_field(line, 47L, 54L, record, line_number),
		type = element,
		index = pdb_integer_field(line, 7L, 11L, record, line_number),
		record = record,
		atom_name = atom_name,
		res_name = res_name,
		chain_id = pdb_character_field(line, 22L, 22L, record, line_number),
		res_seq = pdb_integer_field(line, 23L, 26L, record, line_number),
		i_code = pdb_character_field(line, 27L, 27L, record, line_number),
		alt_loc = pdb_character_field(line, 17L, 17L, record, line_number),
		atom_order = atom_order,
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
parse_pdb_conect_record = function(line, line_number) {
	record = "CONECT"
	line_width = nchar(trimws(line, which = "right"), type = "bytes")
	if (line_width < 11L) {
		stop(
			sprintf(
				"Malformed CONECT record at line %d: expected at least 11 columns",
				line_number
			)
		)
	}

	source = pdb_integer_field(line, 7L, 11L, record, line_number)
	remainder = line_width - 11L

	if (remainder > 0L && remainder %% 5L != 0L) {
		stop(
			sprintf(
				"Malformed CONECT record at line %d: destination fields must be 5 columns wide",
				line_number
			)
		)
	}

	if (remainder == 0L) {
		return(empty_pdb_bonds())
	}

	destinations = c()
	for (start in seq(12L, line_width, by = 5L)) {
		block = substr(line, start, start + 4L)
		if (!nzchar(trimws(block))) {
			next
		}
		destination = suppressWarnings(as.integer(trimws(block)))
		if (is.na(destination)) {
			stop(
				sprintf(
					"Malformed CONECT record at line %d: invalid destination atom serial",
					line_number
				)
			)
		}
		destinations = c(destinations, destination)
	}

	if (length(destinations) == 0L) {
		return(empty_pdb_bonds())
	}

	return(data.frame(
		from = rep(source, length(destinations)),
		to = destinations,
		number = rep(1L, length(destinations)),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
parse_pdb_helix_record = function(line, line_number) {
	record = "HELIX"
	return(data.frame(
		helix_serial = pdb_integer_field(line, 8L, 10L, record, line_number),
		helix_id = trimws(pdb_field(line, 12L, 14L, record, line_number)),
		start_res_name = trimws(pdb_field(line, 16L, 18L, record, line_number)),
		start_chain_id = pdb_character_field(
			line,
			20L,
			20L,
			record,
			line_number
		),
		start_res_seq = pdb_integer_field(line, 22L, 25L, record, line_number),
		start_i_code = pdb_character_field(line, 26L, 26L, record, line_number),
		end_res_name = trimws(pdb_field(line, 28L, 30L, record, line_number)),
		end_chain_id = pdb_character_field(line, 32L, 32L, record, line_number),
		end_res_seq = pdb_integer_field(line, 34L, 37L, record, line_number),
		end_i_code = pdb_character_field(line, 38L, 38L, record, line_number),
		helix_class = pdb_integer_field(line, 39L, 40L, record, line_number),
		comment = trimws(pdb_field(line, 41L, 70L, record, line_number)),
		length = pdb_integer_field(line, 72L, 76L, record, line_number),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
parse_pdb_sheet_record = function(line, line_number) {
	record = "SHEET"
	return(data.frame(
		strand = pdb_integer_field(line, 8L, 10L, record, line_number),
		sheet_id = trimws(pdb_field(line, 12L, 14L, record, line_number)),
		num_strands = pdb_integer_field(line, 15L, 16L, record, line_number),
		start_res_name = trimws(pdb_field(line, 18L, 20L, record, line_number)),
		start_chain_id = pdb_character_field(
			line,
			22L,
			22L,
			record,
			line_number
		),
		start_res_seq = pdb_integer_field(line, 23L, 26L, record, line_number),
		start_i_code = pdb_character_field(line, 27L, 27L, record, line_number),
		end_res_name = trimws(pdb_field(line, 29L, 31L, record, line_number)),
		end_chain_id = pdb_character_field(line, 33L, 33L, record, line_number),
		end_res_seq = pdb_integer_field(line, 34L, 37L, record, line_number),
		end_i_code = pdb_character_field(line, 38L, 38L, record, line_number),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
parse_pdb_ssbond_record = function(line, line_number) {
	record = "SSBOND"
	return(data.frame(
		ser_num = pdb_integer_field(line, 8L, 10L, record, line_number),
		res_name_1 = trimws(pdb_field(line, 12L, 14L, record, line_number)),
		chain_id_1 = pdb_character_field(line, 16L, 16L, record, line_number),
		res_seq_1 = pdb_integer_field(line, 18L, 21L, record, line_number),
		i_code_1 = pdb_character_field(line, 22L, 22L, record, line_number),
		res_name_2 = trimws(pdb_field(line, 26L, 28L, record, line_number)),
		chain_id_2 = pdb_character_field(line, 30L, 30L, record, line_number),
		res_seq_2 = pdb_integer_field(line, 32L, 35L, record, line_number),
		i_code_2 = pdb_character_field(line, 36L, 36L, record, line_number),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
parse_pdb_ter_record = function(line, line_number) {
	record = "TER"
	return(data.frame(
		chain_id = pdb_character_field(line, 22L, 22L, record, line_number),
		res_seq = pdb_integer_field(line, 23L, 26L, record, line_number),
		i_code = pdb_character_field(line, 27L, 27L, record, line_number),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
apply_pdb_secondary_structure = function(residues, helices) {
	if (nrow(helices) == 0L || nrow(residues) == 0L) {
		return(residues)
	}

	for (i in seq_len(nrow(helices))) {
		helix = helices[i, , drop = FALSE]
		if (helix$start_chain_id != helix$end_chain_id) {
			stop("HELIX records spanning multiple chains are not supported yet")
		}
		# Find helix stand and end
		within_range =
			residues$chain_id == helix$start_chain_id &
			pdb_residue_in_range(
				res_seq = residues$res_seq,
				i_code = residues$i_code,
				start_res_seq = helix$start_res_seq,
				start_i_code = helix$start_i_code,
				end_res_seq = helix$end_res_seq,
				end_i_code = helix$end_i_code
			)

		overlapping = within_range & residues$ss_class != "loop"
		if (any(overlapping)) {
			stop(
				"Overlapping HELIX/SHEET residue annotations are not supported yet"
			)
		}

		residues$ss_class[within_range] = "helix"
		residues$helix_id[within_range] = helix$helix_id
	}

	return(residues)
}

#' @keywords internal
apply_pdb_sheet_structure = function(residues, sheets) {
	if (nrow(sheets) == 0L || nrow(residues) == 0L) {
		return(residues)
	}

	for (i in seq_len(nrow(sheets))) {
		sheet = sheets[i, , drop = FALSE]
		if (sheet$start_chain_id != sheet$end_chain_id) {
			stop("SHEET records spanning multiple chains are not supported yet")
		}

		within_range =
			residues$chain_id == sheet$start_chain_id &
			pdb_residue_in_range(
				res_seq = residues$res_seq,
				i_code = residues$i_code,
				start_res_seq = sheet$start_res_seq,
				start_i_code = sheet$start_i_code,
				end_res_seq = sheet$end_res_seq,
				end_i_code = sheet$end_i_code
			)

		overlapping = within_range & residues$ss_class != "loop"
		if (any(overlapping)) {
			stop(
				"Overlapping HELIX/SHEET residue annotations are not supported yet"
			)
		}

		residues$ss_class[within_range] = "sheet"
		residues$sheet_id[within_range] = sheet$sheet_id
	}

	return(residues)
}

#' @keywords internal
select_residue_atom = function(residue_atoms, atom_name) {
	atom_rows = residue_atoms[
		trimws(residue_atoms$atom_name) == atom_name,
		,
		drop = FALSE
	]
	if (nrow(atom_rows) == 0L) {
		return(NULL)
	}

	alt_rank = ifelse(
		atom_rows$alt_loc == "",
		0L,
		ifelse(atom_rows$alt_loc == "A", 1L, 2L)
	)
	atom_rows = atom_rows[order(alt_rank, atom_rows$atom_order), , drop = FALSE]
	return(atom_rows[1, , drop = FALSE])
}

#' @keywords internal
pdb_field = function(line, start, end, record, line_number) {
	if (nchar(line, type = "bytes") < end) {
		stop(
			sprintf(
				"Malformed %s record at line %d: expected at least %d columns",
				record,
				line_number,
				end
			)
		)
	}
	return(substr(line, start, end))
}

#' @keywords internal
pdb_integer_field = function(line, start, end, record, line_number) {
	value = trimws(pdb_field(line, start, end, record, line_number))
	if (!nzchar(value)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: blank integer field in columns %d-%d",
				record,
				line_number,
				start,
				end
			)
		)
	}

	parsed = suppressWarnings(as.integer(value))
	if (is.na(parsed)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: invalid integer field in columns %d-%d",
				record,
				line_number,
				start,
				end
			)
		)
	}
	return(parsed)
}

#' @keywords internal
pdb_numeric_field = function(line, start, end, record, line_number) {
	value = trimws(pdb_field(line, start, end, record, line_number))
	if (!nzchar(value)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: blank numeric field in columns %d-%d",
				record,
				line_number,
				start,
				end
			)
		)
	}

	parsed = suppressWarnings(as.numeric(value))
	if (is.na(parsed)) {
		stop(
			sprintf(
				"Malformed %s record at line %d: invalid numeric field in columns %d-%d",
				record,
				line_number,
				start,
				end
			)
		)
	}
	return(parsed)
}

#' @keywords internal
pdb_character_field = function(line, start, end, record, line_number) {
	return(trimws(pdb_field(line, start, end, record, line_number)))
}

#' @keywords internal
pdb_residue_in_range = function(
	res_seq,
	i_code,
	start_res_seq,
	start_i_code,
	end_res_seq,
	end_i_code
) {
	start_cmp = compare_pdb_residue_key(
		res_seq,
		i_code,
		start_res_seq,
		start_i_code
	)
	end_cmp = compare_pdb_residue_key(res_seq, i_code, end_res_seq, end_i_code)
	return(start_cmp >= 0L & end_cmp <= 0L)
}

#' @keywords internal
compare_pdb_residue_key = function(
	res_seq,
	i_code,
	target_res_seq,
	target_i_code
) {
	result = integer(length(res_seq))
	before = res_seq < target_res_seq
	after = res_seq > target_res_seq
	equal = !before & !after

	result[before] = -1L
	result[after] = 1L

	if (any(equal)) {
		left = normalize_pdb_i_code(i_code[equal])
		right = normalize_pdb_i_code(rep(target_i_code, sum(equal)))
		result[equal] = ifelse(left < right, -1L, ifelse(left > right, 1L, 0L))
	}

	return(result)
}

#' @keywords internal
normalize_pdb_i_code = function(i_code) {
	normalized = ifelse(is.na(i_code) | i_code == "", " ", i_code)
	return(normalized)
}

#' @keywords internal
format_pdb_residue_id = function(chain_id, res_seq, i_code) {
	return(sprintf("%s:%d:%s", chain_id, res_seq, i_code))
}

#' @keywords internal
bind_pdb_rows = function(rows, template) {
	if (length(rows) == 0L) {
		return(template)
	}
	data = do.call(rbind, rows)
	rownames(data) = NULL
	return(data)
}

#' @keywords internal
empty_pdb_atoms = function() {
	return(data.frame(
		x = numeric(),
		y = numeric(),
		z = numeric(),
		type = character(),
		index = integer(),
		record = character(),
		atom_name = character(),
		res_name = character(),
		chain_id = character(),
		res_seq = integer(),
		i_code = character(),
		alt_loc = character(),
		atom_order = integer(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_bonds = function() {
	return(data.frame(
		from = integer(),
		to = integer(),
		number = integer(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_helices = function() {
	return(data.frame(
		helix_serial = integer(),
		helix_id = character(),
		start_res_name = character(),
		start_chain_id = character(),
		start_res_seq = integer(),
		start_i_code = character(),
		end_res_name = character(),
		end_chain_id = character(),
		end_res_seq = integer(),
		end_i_code = character(),
		helix_class = integer(),
		comment = character(),
		length = integer(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_sheets = function() {
	return(data.frame(
		strand = integer(),
		sheet_id = character(),
		num_strands = integer(),
		start_res_name = character(),
		start_chain_id = character(),
		start_res_seq = integer(),
		start_i_code = character(),
		end_res_name = character(),
		end_chain_id = character(),
		end_res_seq = integer(),
		end_i_code = character(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_ssbonds = function() {
	return(data.frame(
		ser_num = integer(),
		res_name_1 = character(),
		chain_id_1 = character(),
		res_seq_1 = integer(),
		i_code_1 = character(),
		res_name_2 = character(),
		chain_id_2 = character(),
		res_seq_2 = integer(),
		i_code_2 = character(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_ters = function() {
	return(data.frame(
		chain_id = character(),
		res_seq = integer(),
		i_code = character(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_residues = function() {
	return(data.frame(
		chain_id = character(),
		res_seq = integer(),
		i_code = character(),
		res_name = character(),
		residue_id = character(),
		n_x = numeric(),
		n_y = numeric(),
		n_z = numeric(),
		ca_x = numeric(),
		ca_y = numeric(),
		ca_z = numeric(),
		c_x = numeric(),
		c_y = numeric(),
		c_z = numeric(),
		o_x = numeric(),
		o_y = numeric(),
		o_z = numeric(),
		has_n = logical(),
		has_ca = logical(),
		has_c = logical(),
		has_o = logical(),
		chain_break_before = logical(),
		chain_break_after = logical(),
		ss_class = character(),
		helix_id = character(),
		sheet_id = character(),
		stringsAsFactors = FALSE
	))
}

#' @keywords internal
empty_pdb_chains = function() {
	return(data.frame(
		chain_id = character(),
		chain_index = integer(),
		start_residue_id = character(),
		end_residue_id = character(),
		terminated_by_ter = logical(),
		stringsAsFactors = FALSE
	))
}

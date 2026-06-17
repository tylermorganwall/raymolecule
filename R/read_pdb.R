#' Read PDB File
#'
#' Reads a legacy fixed-width PDB file and extracts atom locations, bond
#' connections, backbone residue information, and secondary structure records.
#' `model$atoms` and `model$bonds` remain compatible with the existing atom and
#' bond scene generators. Biological assemblies can be generated from
#' `REMARK 350 BIOMT` transforms when requested.
#'
#' @param filename Path to the PDB file.
#' @param atom Default `TRUE`. Whether to include standard residue `ATOM`
#'   records in `model$atoms`.
#' @param nsr Default `TRUE`. Whether to include `HETATM` records in
#'   `model$atoms`.
#' @param assembly Default `"asymmetric_unit"`. Either the deposited asymmetric
#'   unit or the biological assembly generated from `REMARK 350 BIOMT`
#'   transforms.
#' @param assembly_id Default `1L`. Biological assembly identifier to use when
#'   `assembly = "biological"`.
#' @param verbose Default `FALSE`. If `TRUE`, report parsed PDB metadata and
#'   atom/residue/model counts.
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
#'
#'   read_pdb("3nir.pdb", assembly = "biological") |>
#'     generate_ribbon_scene(pathtrace = FALSE) |>
#'     render_model()
#' }
read_pdb = function(
  filename,
  atom = TRUE,
  nsr = TRUE,
  assembly = c("asymmetric_unit", "biological"),
  assembly_id = 1L,
  verbose = FALSE
) {
  assembly = match.arg(assembly)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  if (
    length(assembly_id) != 1L ||
      is.na(assembly_id) ||
      !is.finite(assembly_id) ||
      assembly_id < 1
  ) {
    stop("assembly_id must be a positive integer")
  }
  assembly_id = as.integer(assembly_id)

  #Read in all the lines first
  lines = readLines(filename, warn = FALSE)
  metadata = parse_pdb_metadata(lines, filename)
  #Parse the file
  records = parse_pdb_records(lines)
  if (identical(assembly, "biological")) {
    records = expand_pdb_biological_assembly(
      records = records,
      assembly_id = assembly_id
    )
  }

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
      "serial",
      "record",
      "atom_name",
      "res_name",
      "chain_id",
      "res_seq",
      "i_code",
      "alt_loc",
      "model_id",
      "model_index"
    ),
    drop = FALSE
  ]

  model = list(
    atoms = atoms,
    bonds = records$bonds,
    residues = records$residues,
    helices = records$helices,
    sheets = records$sheets,
    ssbonds = records$ssbonds,
    chains = records$chains,
    assemblies = records$assemblies,
    assembly_mode = assembly,
    assembly_id = if (identical(assembly, "biological")) assembly_id else
      NA_integer_,
    metadata = metadata,
    pdb_type = "pdb"
  )
  if (verbose) {
    message(read_pdb_summary_message(model))
  }

  return(model)
}

#' @keywords internal
parse_pdb_metadata = function(lines, filename = NA_character_) {
  header_line = lines[startsWith(lines, "HEADER")]
  header_line = if (length(header_line) > 0L) header_line[[1]] else ""
  pdb_id = trimws(substr(header_line, 63L, 66L))
  classification = trimws(substr(header_line, 11L, 50L))
  deposition_date = trimws(substr(header_line, 51L, 59L))

  title = parse_pdb_continuation_text(lines, "TITLE", 11L, 80L)
  experiment = parse_pdb_single_record_text(lines, "EXPDTA", 11L, 80L)
  declared_model_text = parse_pdb_single_record_text(lines, "NUMMDL", 11L, 14L)
  declared_model_count = suppressWarnings(as.integer(declared_model_text))

  name = title
  if (!nzchar(name) && nzchar(pdb_id)) {
    name = pdb_id
  }
  if (!nzchar(name) && !is.na(filename)) {
    name = tools::file_path_sans_ext(basename(filename))
  }

  list(
    filename = filename,
    pdb_id = if (nzchar(pdb_id)) pdb_id else NA_character_,
    name = name,
    title = if (nzchar(title)) title else NA_character_,
    classification = if (nzchar(classification)) classification else
      NA_character_,
    deposition_date = if (nzchar(deposition_date)) deposition_date else
      NA_character_,
    experiment = if (nzchar(experiment)) experiment else NA_character_,
    declared_model_count = declared_model_count
  )
}

#' @keywords internal
parse_pdb_single_record_text = function(lines, record, start, end) {
  record_lines = lines[startsWith(lines, record)]
  if (length(record_lines) == 0L) {
    return("")
  }
  trimws(substr(record_lines[[1]], start, end))
}

#' @keywords internal
parse_pdb_continuation_text = function(lines, record, start, end) {
  record_lines = lines[startsWith(lines, record)]
  if (length(record_lines) == 0L) {
    return("")
  }

  title = paste(
    trimws(substr(record_lines, start, end)),
    collapse = " "
  )
  trimws(gsub("[[:space:]]+", " ", title))
}

#' @keywords internal
read_pdb_summary_message = function(model) {
  metadata = model$metadata
  lines = sprintf(
    "Read %s PDB models %s",
    pdb_display_name(model),
    format_pdb_model_set(pdb_model_labels(model))
  )
  if (!is.na(metadata$pdb_id)) {
    lines = c(lines, sprintf("  PDB ID: %s", metadata$pdb_id))
  }
  if (!is.na(metadata$experiment)) {
    lines = c(lines, sprintf("  Experiment: %s", metadata$experiment))
  }
  if (!is.na(metadata$declared_model_count)) {
    lines = c(
      lines,
      sprintf("  Declared models: %d", metadata$declared_model_count)
    )
  }
  lines = c(
    lines,
    sprintf(
      "  Parsed: %d atoms, %d residues, %d chains, %d bonds",
      nrow(model$atoms),
      nrow(model$residues),
      nrow(model$chains),
      nrow(model$bonds)
    )
  )
  if (!identical(model$assembly_mode, "asymmetric_unit")) {
    lines = c(
      lines,
      sprintf(
        "  Assembly: %s %d",
        model$assembly_mode,
        model$assembly_id
      )
    )
  }

  paste(lines, collapse = "\n")
}

#' @keywords internal
pdb_display_name = function(model) {
  if (
    !is.null(model$metadata) &&
      !is.na(model$metadata$name) &&
      nzchar(model$metadata$name)
  ) {
    return(model$metadata$name)
  }
  if (
    !is.null(model$metadata) &&
      !is.na(model$metadata$pdb_id) &&
      nzchar(model$metadata$pdb_id)
  ) {
    return(model$metadata$pdb_id)
  }
  "unknown"
}

#' @keywords internal
pdb_model_labels = function(model) {
  source = NULL
  if (!is.null(model$residues) && nrow(model$residues) > 0L) {
    source = model$residues
  } else if (!is.null(model$atoms) && nrow(model$atoms) > 0L) {
    source = model$atoms
  }
  if (is.null(source) || !"model_index" %in% names(source)) {
    return(NA_integer_)
  }

  source = source[order(source$model_index), , drop = FALSE]
  model_rows = !duplicated(source$model_index)
  model_ids = source$model_id[model_rows]
  model_indices = source$model_index[model_rows]
  if (all(is.na(model_ids))) {
    return(model_indices)
  }
  ifelse(is.na(model_ids), model_indices, model_ids)
}

#' @keywords internal
format_pdb_model_set = function(model_ids) {
  model_ids = as.integer(model_ids)
  model_ids = model_ids[!is.na(model_ids)]
  if (length(model_ids) == 0L) {
    return("[unknown]")
  }
  model_ids = sort(unique(model_ids))
  runs = list()
  run_start = model_ids[[1]]
  previous = model_ids[[1]]
  counter = 1L

  if (length(model_ids) > 1L) {
    for (model_id in model_ids[-1]) {
      if (model_id == previous + 1L) {
        previous = model_id
        next
      }
      runs[[counter]] = c(run_start, previous)
      counter = counter + 1L
      run_start = model_id
      previous = model_id
    }
  }
  runs[[counter]] = c(run_start, previous)

  paste0(
    "[",
    paste(
      vapply(
        runs,
        function(run) {
          if (identical(run[[1]], run[[2]])) {
            return(as.character(run[[1]]))
          }
          sprintf("%d-%d", run[[1]], run[[2]])
        },
        character(1)
      ),
      collapse = ", "
    ),
    "]"
  )
}

#' @keywords internal
parse_pdb_records = function(lines) {
  #Set up lists for all different objects within the ribbon to render
  atom_lines = list()
  atom_line_numbers = integer()
  atom_records = character()
  atom_orders = integer()
  atom_model_ids = integer()
  atom_model_indices = integer()
  bond_rows = list()
  helix_rows = list()
  sheet_rows = list()
  ssbond_rows = list()
  ter_rows = list()
  assembly_rows = list()

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
  current_model_id = NA_integer_
  current_model_index = 1L

  # Parse each line
  for (line_number in seq_along(lines)) {
    line = lines[[line_number]]
    #Load the first 6 characters to get the record type
    record = trimws(substr(line, 1L, 6L))

    if (identical(record, "MODEL")) {
      model_count = model_count + 1L
      current_model_id = parse_pdb_model_id(
        line = line,
        model_index = model_count
      )
      current_model_index = model_count
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
      atom_lines[[atom_counter]] = line
      atom_line_numbers[[atom_counter]] = line_number
      atom_records[[atom_counter]] = record
      atom_orders[[atom_counter]] = atom_order
      atom_model_ids[[atom_counter]] = current_model_id
      atom_model_indices[[atom_counter]] = current_model_index
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
        line_number = line_number,
        model_id = current_model_id,
        model_index = current_model_index
      )
      ter_counter = ter_counter + 1L
      next
    }
  }

  assembly_rows = parse_pdb_assembly_records(lines)

  # This specifies the schema for each data type, and ensures it's enforced even at zero rows
  atoms = parse_pdb_atom_records(
    lines = unlist(atom_lines, use.names = FALSE),
    line_numbers = atom_line_numbers,
    records = atom_records,
    atom_orders = atom_orders,
    model_ids = atom_model_ids,
    model_indices = atom_model_indices
  )
  bonds = bind_pdb_rows(bond_rows, empty_pdb_bonds())
  normalized_indices = normalize_pdb_atom_and_bond_indices(
    atoms = atoms,
    bonds = bonds
  )
  atoms = normalized_indices$atoms
  bonds = normalized_indices$bonds
  helices = bind_pdb_rows(helix_rows, empty_pdb_helices())
  sheets = bind_pdb_rows(sheet_rows, empty_pdb_sheets())
  ssbonds = bind_pdb_rows(ssbond_rows, empty_pdb_ssbonds())
  ters = bind_pdb_rows(ter_rows, empty_pdb_ters())
  assemblies = bind_pdb_rows(assembly_rows, empty_pdb_assemblies())

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
    ssbonds = ssbonds,
    ters = ters,
    assemblies = assemblies
  ))
}

#' @keywords internal
parse_pdb_model_id = function(line, model_index) {
  value = trimws(substr(line, 11L, 14L))
  if (!nzchar(value)) {
    return(as.integer(model_index))
  }

  parsed = suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    return(as.integer(model_index))
  }
  return(parsed)
}

#' @keywords internal
normalize_pdb_atom_and_bond_indices = function(atoms, bonds) {
  if (nrow(atoms) == 0L) {
    bonds = empty_pdb_bonds()
    return(list(atoms = atoms, bonds = bonds))
  }

  if (anyDuplicated(atoms$index) > 0L) {
    atoms$index = seq_len(nrow(atoms))
  }

  if (nrow(bonds) == 0L) {
    bonds = empty_pdb_bonds()
    return(list(atoms = atoms, bonds = bonds))
  }

  model_indices = unique(atoms$model_index)
  bond_rows = list()
  bond_counter = 1L

  for (model_index in model_indices) {
    model_atoms = atoms[atoms$model_index == model_index, , drop = FALSE]
    if (nrow(model_atoms) == 0L) {
      next
    }

    serial_keys = as.character(model_atoms$serial)
    unique_serial = !duplicated(serial_keys)
    atom_map = stats::setNames(
      model_atoms$index[unique_serial],
      serial_keys[unique_serial]
    )
    selected_bonds = bonds[
      as.character(bonds$from) %in%
        names(atom_map) &
        as.character(bonds$to) %in% names(atom_map),
      ,
      drop = FALSE
    ]
    if (nrow(selected_bonds) == 0L) {
      next
    }

    selected_bonds$from = unname(atom_map[as.character(selected_bonds$from)])
    selected_bonds$to = unname(atom_map[as.character(selected_bonds$to)])
    selected_bonds$model_id = model_atoms$model_id[1]
    selected_bonds$model_index = model_index
    bond_rows[[bond_counter]] = selected_bonds
    bond_counter = bond_counter + 1L
  }

  bonds = bind_pdb_rows(bond_rows, empty_pdb_bonds())
  return(list(atoms = atoms, bonds = bonds))
}

#' @keywords internal
parse_pdb_assembly_records = function(lines) {
  assembly_rows = list()
  current_assembly_ids = integer()
  current_chain_ids = character()
  current_chain_block_id = 0L
  matrix_registry = list()
  row_counter = 1L

  for (line_number in seq_along(lines)) {
    line = lines[[line_number]]
    if (!startsWith(line, "REMARK 350")) {
      next
    }

    if (grepl("^REMARK 350 BIOMOLECULE:", line)) {
      current_assembly_ids = parse_pdb_biomolecule_ids(line, line_number)
      current_chain_ids = character()
      current_chain_block_id = 0L
      next
    }

    if (grepl("^REMARK 350 APPLY THE FOLLOWING TO CHAINS:", line)) {
      current_chain_ids = parse_pdb_assembly_chain_ids(line)
      current_chain_block_id = current_chain_block_id + 1L
      next
    }

    if (grepl("^REMARK 350                    AND CHAINS:", line)) {
      current_chain_ids = c(
        current_chain_ids,
        parse_pdb_assembly_chain_ids(line)
      )
      current_chain_ids = unique(current_chain_ids)
      next
    }

    if (grepl("^REMARK 350   BIOMT[123]", line)) {
      if (
        length(current_assembly_ids) == 0L || length(current_chain_ids) == 0L
      ) {
        stop(
          sprintf(
            "Malformed REMARK 350 BIOMT record at line %d: missing BIOMOLECULE or chain assignment",
            line_number
          )
        )
      }

      biomt_row = parse_pdb_biomt_record(line, line_number)
      for (assembly_id in current_assembly_ids) {
        matrix_key = paste(
          assembly_id,
          current_chain_block_id,
          biomt_row$transform_id,
          sep = "\r"
        )
        if (is.null(matrix_registry[[matrix_key]])) {
          matrix_registry[[matrix_key]] = list(
            rows = vector(mode = "list", length = 3L),
            chains = current_chain_ids
          )
        }

        matrix_registry[[matrix_key]]$rows[[
          biomt_row$row_index
        ]] = biomt_row$values
        matrix_registry[[matrix_key]]$chains = current_chain_ids

        if (
          all(vapply(
            matrix_registry[[matrix_key]]$rows,
            Negate(is.null),
            logical(1)
          ))
        ) {
          matrix_values = do.call(rbind, matrix_registry[[matrix_key]]$rows)
          for (chain_id in matrix_registry[[matrix_key]]$chains) {
            assembly_rows[[row_counter]] = data.frame(
              assembly_id = assembly_id,
              transform_id = biomt_row$transform_id,
              chain_id = chain_id,
              m11 = matrix_values[1, 1],
              m12 = matrix_values[1, 2],
              m13 = matrix_values[1, 3],
              t1 = matrix_values[1, 4],
              m21 = matrix_values[2, 1],
              m22 = matrix_values[2, 2],
              m23 = matrix_values[2, 3],
              t2 = matrix_values[2, 4],
              m31 = matrix_values[3, 1],
              m32 = matrix_values[3, 2],
              m33 = matrix_values[3, 3],
              t3 = matrix_values[3, 4],
              stringsAsFactors = FALSE
            )
            row_counter = row_counter + 1L
          }
        }
      }
    }
  }

  if (length(matrix_registry) > 0L) {
    complete_flags = vapply(
      matrix_registry,
      function(entry) all(vapply(entry$rows, Negate(is.null), logical(1))),
      logical(1)
    )
    if (any(!complete_flags)) {
      stop("Malformed REMARK 350 BIOMT records: incomplete transform matrix")
    }
  }

  return(assembly_rows)
}

#' @keywords internal
parse_pdb_biomolecule_ids = function(line, line_number) {
  text = trimws(sub("^REMARK 350 BIOMOLECULE:\\s*", "", line))
  if (!nzchar(text)) {
    stop(
      sprintf(
        "Malformed REMARK 350 BIOMOLECULE record at line %d",
        line_number
      )
    )
  }

  ids = as.integer(trimws(unlist(strsplit(text, ","))))
  if (any(is.na(ids))) {
    stop(
      sprintf(
        "Malformed REMARK 350 BIOMOLECULE record at line %d",
        line_number
      )
    )
  }

  return(ids)
}

#' @keywords internal
parse_pdb_assembly_chain_ids = function(line) {
  text = trimws(sub(
    "^REMARK 350\\s+(APPLY THE FOLLOWING TO CHAINS:|AND CHAINS:)\\s*",
    "",
    line
  ))
  chain_ids = trimws(unlist(strsplit(text, ",")))
  chain_ids = chain_ids[nzchar(chain_ids)]
  return(chain_ids)
}

#' @keywords internal
parse_pdb_biomt_record = function(line, line_number) {
  matches = regexec(
    "^REMARK 350   BIOMT([123])\\s+(\\d+)\\s+([-0-9.Ee+]+)\\s+([-0-9.Ee+]+)\\s+([-0-9.Ee+]+)\\s+([-0-9.Ee+]+)\\s*$",
    line
  )
  values = regmatches(line, matches)[[1]]
  if (length(values) != 7L) {
    stop(
      sprintf(
        "Malformed REMARK 350 BIOMT record at line %d",
        line_number
      )
    )
  }

  parsed_values = as.numeric(values[4:7])
  if (any(!is.finite(parsed_values))) {
    stop(
      sprintf(
        "Malformed REMARK 350 BIOMT record at line %d",
        line_number
      )
    )
  }

  return(list(
    row_index = as.integer(values[2]),
    transform_id = as.integer(values[3]),
    values = parsed_values
  ))
}

#' @keywords internal
expand_pdb_biological_assembly = function(records, assembly_id) {
  if (!"assemblies" %in% names(records) || nrow(records$assemblies) == 0L) {
    stop("Biological assembly transforms were not found in this PDB file")
  }

  assembly_rows = records$assemblies[
    records$assemblies$assembly_id == assembly_id,
    ,
    drop = FALSE
  ]
  if (nrow(assembly_rows) == 0L) {
    stop(sprintf(
      "Biological assembly %d was not found in this PDB file",
      assembly_id
    ))
  }

  transform_ids = unique(assembly_rows$transform_id)
  transform_ids = transform_ids[order(transform_ids)]
  expanded_atom_rows = list()
  expanded_bond_rows = list()
  expanded_helix_rows = list()
  expanded_sheet_rows = list()
  expanded_ssbond_rows = list()
  expanded_ter_rows = list()
  atom_counter = 1L
  bond_counter = 1L
  helix_counter = 1L
  sheet_counter = 1L
  ssbond_counter = 1L
  ter_counter = 1L
  next_atom_index = 1L
  next_atom_order = 1L

  for (transform_position in seq_along(transform_ids)) {
    transform_id = transform_ids[[transform_position]]
    transform_rows = assembly_rows[
      assembly_rows$transform_id == transform_id,
      ,
      drop = FALSE
    ]
    source_chain_ids = transform_rows$chain_id
    source_atoms = records$atoms[
      records$atoms$chain_id %in% source_chain_ids,
      ,
      drop = FALSE
    ]
    if (nrow(source_atoms) == 0L) {
      next
    }

    transform_map = integer()
    for (row_index in seq_len(nrow(transform_rows))) {
      transform_row = transform_rows[row_index, , drop = FALSE]
      source_chain_id = transform_row$chain_id
      target_chain_id = format_pdb_assembly_chain_id(
        source_chain_id = source_chain_id,
        transform_position = transform_position
      )
      chain_atoms = records$atoms[
        records$atoms$chain_id == source_chain_id,
        ,
        drop = FALSE
      ]
      if (nrow(chain_atoms) == 0L) {
        next
      }

      rotation = matrix(
        c(
          transform_row$m11,
          transform_row$m12,
          transform_row$m13,
          transform_row$m21,
          transform_row$m22,
          transform_row$m23,
          transform_row$m31,
          transform_row$m32,
          transform_row$m33
        ),
        nrow = 3L,
        byrow = TRUE
      )
      translation = c(transform_row$t1, transform_row$t2, transform_row$t3)
      coords = as.matrix(chain_atoms[, c("x", "y", "z"), drop = FALSE])
      transformed_coords = coords %*%
        t(rotation) +
        matrix(rep(translation, each = nrow(coords)), ncol = 3L)

      transformed_atoms = chain_atoms
      transformed_atoms$x = transformed_coords[, 1]
      transformed_atoms$y = transformed_coords[, 2]
      transformed_atoms$z = transformed_coords[, 3]
      transformed_atoms$chain_id = target_chain_id
      transformed_atoms$index = seq.int(
        from = next_atom_index,
        length.out = nrow(transformed_atoms)
      )
      transformed_atoms$atom_order = seq.int(
        from = next_atom_order,
        length.out = nrow(transformed_atoms)
      )

      transform_map = c(
        transform_map,
        stats::setNames(transformed_atoms$index, chain_atoms$index)
      )
      next_atom_index = next_atom_index + nrow(transformed_atoms)
      next_atom_order = next_atom_order + nrow(transformed_atoms)
      expanded_atom_rows[[atom_counter]] = transformed_atoms
      atom_counter = atom_counter + 1L

      chain_helices = records$helices[
        records$helices$start_chain_id == source_chain_id &
          records$helices$end_chain_id == source_chain_id,
        ,
        drop = FALSE
      ]
      if (nrow(chain_helices) > 0L) {
        chain_helices$start_chain_id = target_chain_id
        chain_helices$end_chain_id = target_chain_id
        expanded_helix_rows[[helix_counter]] = chain_helices
        helix_counter = helix_counter + 1L
      }

      chain_sheets = records$sheets[
        records$sheets$start_chain_id == source_chain_id &
          records$sheets$end_chain_id == source_chain_id,
        ,
        drop = FALSE
      ]
      if (nrow(chain_sheets) > 0L) {
        chain_sheets$start_chain_id = target_chain_id
        chain_sheets$end_chain_id = target_chain_id
        expanded_sheet_rows[[sheet_counter]] = chain_sheets
        sheet_counter = sheet_counter + 1L
      }

      chain_ters = records$ters[
        records$ters$chain_id == source_chain_id,
        ,
        drop = FALSE
      ]
      if (nrow(chain_ters) > 0L) {
        chain_ters$chain_id = target_chain_id
        expanded_ter_rows[[ter_counter]] = chain_ters
        ter_counter = ter_counter + 1L
      }
    }

    if (length(transform_map) > 0L) {
      selected_bonds = records$bonds[
        as.character(records$bonds$from) %in%
          names(transform_map) &
          as.character(records$bonds$to) %in% names(transform_map),
        ,
        drop = FALSE
      ]
      if (nrow(selected_bonds) > 0L) {
        selected_bonds$from = unname(transform_map[as.character(
          selected_bonds$from
        )])
        selected_bonds$to = unname(transform_map[as.character(
          selected_bonds$to
        )])
        expanded_bond_rows[[bond_counter]] = selected_bonds
        bond_counter = bond_counter + 1L
      }
    }

    transform_ssbonds = records$ssbonds[
      records$ssbonds$chain_id_1 %in%
        source_chain_ids &
        records$ssbonds$chain_id_2 %in% source_chain_ids,
      ,
      drop = FALSE
    ]
    if (nrow(transform_ssbonds) > 0L) {
      transform_ssbonds$chain_id_1 = vapply(
        transform_ssbonds$chain_id_1,
        format_pdb_assembly_chain_id,
        character(1),
        transform_position = transform_position
      )
      transform_ssbonds$chain_id_2 = vapply(
        transform_ssbonds$chain_id_2,
        format_pdb_assembly_chain_id,
        character(1),
        transform_position = transform_position
      )
      expanded_ssbond_rows[[ssbond_counter]] = transform_ssbonds
      ssbond_counter = ssbond_counter + 1L
    }
  }

  expanded_atoms = bind_pdb_rows(expanded_atom_rows, empty_pdb_atoms())
  expanded_bonds = bind_pdb_rows(expanded_bond_rows, empty_pdb_bonds())
  expanded_helices = bind_pdb_rows(expanded_helix_rows, empty_pdb_helices())
  expanded_sheets = bind_pdb_rows(expanded_sheet_rows, empty_pdb_sheets())
  expanded_ssbonds = bind_pdb_rows(expanded_ssbond_rows, empty_pdb_ssbonds())
  expanded_ters = bind_pdb_rows(expanded_ter_rows, empty_pdb_ters())
  expanded_residue_info = build_pdb_residues(
    atoms = expanded_atoms[expanded_atoms$record == "ATOM", , drop = FALSE],
    helices = expanded_helices,
    sheets = expanded_sheets,
    ters = expanded_ters
  )

  return(list(
    atoms = expanded_atoms,
    bonds = expanded_bonds,
    residues = expanded_residue_info$residues,
    chains = expanded_residue_info$chains,
    helices = expanded_helices,
    sheets = expanded_sheets,
    ssbonds = expanded_ssbonds,
    ters = expanded_ters,
    assemblies = records$assemblies
  ))
}

#' @keywords internal
format_pdb_assembly_chain_id = function(source_chain_id, transform_position) {
  if (transform_position == 1L) {
    return(source_chain_id)
  }
  return(sprintf("%s_%d", source_chain_id, transform_position))
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
    atoms$model_index,
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
      model_id = first_atom$model_id,
      model_index = first_atom$model_index,
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
      sheet_strand = NA_integer_,
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
    chain_change =
      residues$model_index[-1] != residues$model_index[-nrow(residues)] |
      residues$chain_id[-1] != residues$chain_id[-nrow(residues)]
    # Mark the ends
    residue_break_after[-nrow(residues)] =
      residue_break_after[-nrow(residues)] | chain_change
    # Mark the beginnings
    residue_break_before[-1] =
      residue_break_before[-1] | chain_change
  }

  if (nrow(ters) > 0L) {
    residue_match = match(
      paste(
        ters$model_index,
        ters$chain_id,
        ters$res_seq,
        ters$i_code,
        sep = "\r"
      ),
      paste(
        residues$model_index,
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
      model_id = residues$model_id[chain_start],
      model_index = residues$model_index[chain_start],
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
parse_pdb_atom_records = function(
  lines,
  line_numbers,
  records,
  atom_orders,
  model_ids,
  model_indices
) {
  if (length(lines) == 0L) {
    return(empty_pdb_atoms())
  }

  atom_name = trimws(pdb_vector_field(lines, 13L, 16L, records, line_numbers))
  blank_atom_name = which(!nzchar(atom_name))
  if (length(blank_atom_name) > 0L) {
    i = blank_atom_name[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: blank atom name in columns 13-16",
      records[[i]],
      line_numbers[[i]]
    ))
  }

  res_name = trimws(pdb_vector_field(lines, 18L, 20L, records, line_numbers))
  blank_res_name = which(!nzchar(res_name))
  if (length(blank_res_name) > 0L) {
    i = blank_res_name[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: blank residue name in columns 18-20",
      records[[i]],
      line_numbers[[i]]
    ))
  }

  element = trimws(pdb_vector_field(lines, 77L, 78L, records, line_numbers))
  blank_element = which(!nzchar(element))
  if (length(blank_element) > 0L) {
    i = blank_element[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: blank element field in columns 77-78",
      records[[i]],
      line_numbers[[i]]
    ))
  }

  serial = pdb_integer_values(lines, 7L, 11L, records, line_numbers)

  return(data.frame(
    x = pdb_numeric_values(lines, 31L, 38L, records, line_numbers),
    y = pdb_numeric_values(lines, 39L, 46L, records, line_numbers),
    z = pdb_numeric_values(lines, 47L, 54L, records, line_numbers),
    type = element,
    index = serial,
    serial = serial,
    record = records,
    atom_name = atom_name,
    res_name = res_name,
    chain_id = trimws(pdb_vector_field(lines, 22L, 22L, records, line_numbers)),
    res_seq = pdb_integer_values(lines, 23L, 26L, records, line_numbers),
    i_code = trimws(pdb_vector_field(lines, 27L, 27L, records, line_numbers)),
    alt_loc = trimws(pdb_vector_field(lines, 17L, 17L, records, line_numbers)),
    model_id = model_ids,
    model_index = model_indices,
    atom_order = atom_orders,
    stringsAsFactors = FALSE
  ))
}

#' @keywords internal
parse_pdb_atom_record = function(
  line,
  line_number,
  record,
  atom_order,
  model_id,
  model_index
) {
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

  serial = pdb_integer_field(line, 7L, 11L, record, line_number)

  return(data.frame(
    x = pdb_numeric_field(line, 31L, 38L, record, line_number),
    y = pdb_numeric_field(line, 39L, 46L, record, line_number),
    z = pdb_numeric_field(line, 47L, 54L, record, line_number),
    type = element,
    index = serial,
    serial = serial,
    record = record,
    atom_name = atom_name,
    res_name = res_name,
    chain_id = pdb_character_field(line, 22L, 22L, record, line_number),
    res_seq = pdb_integer_field(line, 23L, 26L, record, line_number),
    i_code = pdb_character_field(line, 27L, 27L, record, line_number),
    alt_loc = pdb_character_field(line, 17L, 17L, record, line_number),
    model_id = model_id,
    model_index = model_index,
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
parse_pdb_ter_record = function(line, line_number, model_id, model_index) {
  record = "TER"
  return(data.frame(
    chain_id = pdb_character_field(line, 22L, 22L, record, line_number),
    res_seq = pdb_integer_field(line, 23L, 26L, record, line_number),
    i_code = pdb_character_field(line, 27L, 27L, record, line_number),
    model_id = model_id,
    model_index = model_index,
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

    # Real PDB files can contain overlapping secondary-structure annotations.
    # Keep the existing assignment and only fill unassigned loop residues.
    new_rows = within_range & residues$ss_class == "loop"
    residues$ss_class[new_rows] = "helix"
    residues$helix_id[new_rows] = helix$helix_id
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

    # Real PDB files can contain overlapping secondary-structure annotations.
    # Keep the existing assignment and only fill unassigned loop residues.
    new_rows = within_range & residues$ss_class == "loop"
    residues$ss_class[new_rows] = "sheet"
    residues$sheet_id[new_rows] = sheet$sheet_id
    residues$sheet_strand[new_rows] = sheet$strand
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
pdb_vector_field = function(lines, start, end, records, line_numbers) {
  line_widths = nchar(lines, type = "bytes")
  malformed = which(line_widths < end)
  if (length(malformed) > 0L) {
    i = malformed[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: expected at least %d columns",
      records[[i]],
      line_numbers[[i]],
      end
    ))
  }
  return(substr(lines, start, end))
}

#' @keywords internal
pdb_integer_values = function(lines, start, end, records, line_numbers) {
  value = trimws(pdb_vector_field(lines, start, end, records, line_numbers))
  blank = which(!nzchar(value))
  if (length(blank) > 0L) {
    i = blank[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: blank integer field in columns %d-%d",
      records[[i]],
      line_numbers[[i]],
      start,
      end
    ))
  }

  parsed = suppressWarnings(as.integer(value))
  invalid = which(is.na(parsed))
  if (length(invalid) > 0L) {
    i = invalid[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: invalid integer field in columns %d-%d",
      records[[i]],
      line_numbers[[i]],
      start,
      end
    ))
  }
  return(parsed)
}

#' @keywords internal
pdb_numeric_values = function(lines, start, end, records, line_numbers) {
  value = trimws(pdb_vector_field(lines, start, end, records, line_numbers))
  blank = which(!nzchar(value))
  if (length(blank) > 0L) {
    i = blank[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: blank numeric field in columns %d-%d",
      records[[i]],
      line_numbers[[i]],
      start,
      end
    ))
  }

  parsed = suppressWarnings(as.numeric(value))
  invalid = which(is.na(parsed))
  if (length(invalid) > 0L) {
    i = invalid[[1]]
    stop(sprintf(
      "Malformed %s record at line %d: invalid numeric field in columns %d-%d",
      records[[i]],
      line_numbers[[i]],
      start,
      end
    ))
  }
  return(parsed)
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
    serial = integer(),
    record = character(),
    atom_name = character(),
    res_name = character(),
    chain_id = character(),
    res_seq = integer(),
    i_code = character(),
    alt_loc = character(),
    model_id = integer(),
    model_index = integer(),
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
    model_id = integer(),
    model_index = integer(),
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
    model_id = integer(),
    model_index = integer(),
    stringsAsFactors = FALSE
  ))
}

#' @keywords internal
empty_pdb_assemblies = function() {
  return(data.frame(
    assembly_id = integer(),
    transform_id = integer(),
    chain_id = character(),
    m11 = numeric(),
    m12 = numeric(),
    m13 = numeric(),
    t1 = numeric(),
    m21 = numeric(),
    m22 = numeric(),
    m23 = numeric(),
    t2 = numeric(),
    m31 = numeric(),
    m32 = numeric(),
    m33 = numeric(),
    t3 = numeric(),
    stringsAsFactors = FALSE
  ))
}

#' @keywords internal
empty_pdb_residues = function() {
  return(data.frame(
    chain_id = character(),
    model_id = integer(),
    model_index = integer(),
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
    sheet_strand = integer(),
    stringsAsFactors = FALSE
  ))
}

#' @keywords internal
empty_pdb_chains = function() {
  return(data.frame(
    chain_id = character(),
    model_id = integer(),
    model_index = integer(),
    chain_index = integer(),
    start_residue_id = character(),
    end_residue_id = character(),
    terminated_by_ter = logical(),
    stringsAsFactors = FALSE
  ))
}

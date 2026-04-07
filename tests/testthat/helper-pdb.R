make_pdb_line = function(width = 80) {
  return(paste(rep(" ", width), collapse = ""))
}

set_pdb_field = function(line, start, end, value, align = c("left", "right")) {
  align = match.arg(align)
  width = end - start + 1L
  text = as.character(value)
  formatted = if (align == "right") {
    sprintf(paste0("%", width, "s"), text)
  } else {
    sprintf(paste0("%-", width, "s"), text)
  }
  substring(line, start, end) = substr(formatted, 1L, width)
  return(line)
}

format_atom_line = function(
  serial,
  atom_name,
  res_name,
  chain_id,
  res_seq,
  x,
  y,
  z,
  element,
  record = "ATOM",
  i_code = "",
  alt_loc = ""
) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, record, align = "left")
  line = set_pdb_field(line, 7, 11, serial, align = "right")
  line = set_pdb_field(line, 13, 16, atom_name, align = "left")
  line = set_pdb_field(line, 17, 17, alt_loc, align = "left")
  line = set_pdb_field(line, 18, 20, res_name, align = "left")
  line = set_pdb_field(line, 22, 22, chain_id, align = "left")
  line = set_pdb_field(line, 23, 26, res_seq, align = "right")
  line = set_pdb_field(line, 27, 27, i_code, align = "left")
  line = set_pdb_field(line, 31, 38, sprintf("%.3f", x), align = "right")
  line = set_pdb_field(line, 39, 46, sprintf("%.3f", y), align = "right")
  line = set_pdb_field(line, 47, 54, sprintf("%.3f", z), align = "right")
  line = set_pdb_field(line, 55, 60, sprintf("%.2f", 1), align = "right")
  line = set_pdb_field(line, 61, 66, sprintf("%.2f", 20), align = "right")
  line = set_pdb_field(line, 77, 78, element, align = "right")
  return(line)
}

format_ter_line = function(serial, res_name, chain_id, res_seq, i_code = "") {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "TER", align = "left")
  line = set_pdb_field(line, 7, 11, serial, align = "right")
  line = set_pdb_field(line, 18, 20, res_name, align = "left")
  line = set_pdb_field(line, 22, 22, chain_id, align = "left")
  line = set_pdb_field(line, 23, 26, res_seq, align = "right")
  line = set_pdb_field(line, 27, 27, i_code, align = "left")
  return(line)
}

format_helix_line = function(
  serial,
  helix_id,
  start_res_name,
  start_chain_id,
  start_res_seq,
  end_res_name,
  end_chain_id,
  end_res_seq,
  start_i_code = "",
  end_i_code = "",
  helix_class = 1,
  length = 2
) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "HELIX", align = "left")
  line = set_pdb_field(line, 8, 10, serial, align = "right")
  line = set_pdb_field(line, 12, 14, helix_id, align = "left")
  line = set_pdb_field(line, 16, 18, start_res_name, align = "left")
  line = set_pdb_field(line, 20, 20, start_chain_id, align = "left")
  line = set_pdb_field(line, 22, 25, start_res_seq, align = "right")
  line = set_pdb_field(line, 26, 26, start_i_code, align = "left")
  line = set_pdb_field(line, 28, 30, end_res_name, align = "left")
  line = set_pdb_field(line, 32, 32, end_chain_id, align = "left")
  line = set_pdb_field(line, 34, 37, end_res_seq, align = "right")
  line = set_pdb_field(line, 38, 38, end_i_code, align = "left")
  line = set_pdb_field(line, 39, 40, helix_class, align = "right")
  line = set_pdb_field(line, 72, 76, length, align = "right")
  return(line)
}

format_sheet_line = function(
  strand,
  sheet_id,
  num_strands,
  start_res_name,
  start_chain_id,
  start_res_seq,
  end_res_name,
  end_chain_id,
  end_res_seq,
  start_i_code = "",
  end_i_code = ""
) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "SHEET", align = "left")
  line = set_pdb_field(line, 8, 10, strand, align = "right")
  line = set_pdb_field(line, 12, 14, sheet_id, align = "left")
  line = set_pdb_field(line, 15, 16, num_strands, align = "right")
  line = set_pdb_field(line, 18, 20, start_res_name, align = "left")
  line = set_pdb_field(line, 22, 22, start_chain_id, align = "left")
  line = set_pdb_field(line, 23, 26, start_res_seq, align = "right")
  line = set_pdb_field(line, 27, 27, start_i_code, align = "left")
  line = set_pdb_field(line, 29, 31, end_res_name, align = "left")
  line = set_pdb_field(line, 33, 33, end_chain_id, align = "left")
  line = set_pdb_field(line, 34, 37, end_res_seq, align = "right")
  line = set_pdb_field(line, 38, 38, end_i_code, align = "left")
  return(line)
}

format_ssbond_line = function(
  serial,
  chain_id_1,
  res_seq_1,
  chain_id_2,
  res_seq_2,
  i_code_1 = "",
  i_code_2 = ""
) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "SSBOND", align = "left")
  line = set_pdb_field(line, 8, 10, serial, align = "right")
  line = set_pdb_field(line, 12, 14, "CYS", align = "left")
  line = set_pdb_field(line, 16, 16, chain_id_1, align = "left")
  line = set_pdb_field(line, 18, 21, res_seq_1, align = "right")
  line = set_pdb_field(line, 22, 22, i_code_1, align = "left")
  line = set_pdb_field(line, 26, 28, "CYS", align = "left")
  line = set_pdb_field(line, 30, 30, chain_id_2, align = "left")
  line = set_pdb_field(line, 32, 35, res_seq_2, align = "right")
  line = set_pdb_field(line, 36, 36, i_code_2, align = "left")
  return(line)
}

format_conect_line = function(source, destinations) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "CONECT", align = "left")
  line = set_pdb_field(line, 7, 11, source, align = "right")
  starts = seq(12, by = 5, length.out = length(destinations))
  for (i in seq_along(destinations)) {
    line = set_pdb_field(line, starts[i], starts[i] + 4L, destinations[[i]], align = "right")
  }
  return(line)
}

format_model_line = function(serial) {
  line = make_pdb_line()
  line = set_pdb_field(line, 1, 6, "MODEL", align = "left")
  line = set_pdb_field(line, 11, 14, serial, align = "right")
  return(line)
}

format_remark350_biomolecule_line = function(assembly_id) {
  return(sprintf("REMARK 350 BIOMOLECULE: %s", assembly_id))
}

format_remark350_apply_chains_line = function(chain_ids) {
  return(sprintf(
    "REMARK 350 APPLY THE FOLLOWING TO CHAINS: %s",
    paste(chain_ids, collapse = ", ")
  ))
}

format_remark350_biomt_line = function(row_index, transform_id, values) {
  if (length(values) != 4L) {
    stop("format_remark350_biomt_line() requires four numeric values")
  }
  return(sprintf(
    "REMARK 350   BIOMT%d %3d %10.6f %10.6f %10.6f %14.5f",
    row_index,
    transform_id,
    values[[1]],
    values[[2]],
    values[[3]],
    values[[4]]
  ))
}

backbone_residue_lines = function(serial_start, chain_id, res_seq, res_name, ca, i_code = "") {
  n = ca + c(-0.6, 0.0, 0.0)
  c_atom = ca + c(0.6, 0.0, 0.0)
  o = c_atom + c(0.2, 0.8, 0.1)

  lines = c(
    format_atom_line(serial_start, "N", res_name, chain_id, res_seq, n[1], n[2], n[3], "N", i_code = i_code),
    format_atom_line(serial_start + 1L, "CA", res_name, chain_id, res_seq, ca[1], ca[2], ca[3], "C", i_code = i_code),
    format_atom_line(serial_start + 2L, "C", res_name, chain_id, res_seq, c_atom[1], c_atom[2], c_atom[3], "C", i_code = i_code),
    format_atom_line(serial_start + 3L, "O", res_name, chain_id, res_seq, o[1], o[2], o[3], "O", i_code = i_code)
  )

  return(list(lines = lines, next_serial = serial_start + 4L))
}

write_pdb_fixture = function(lines) {
  file = tempfile(fileext = ".pdb")
  writeLines(lines, file)
  return(file)
}

mesh_edge_table = function(indices) {
  edges = rbind(indices[, c(1, 2)], indices[, c(2, 3)], indices[, c(3, 1)])
  edges = t(apply(edges, 1, sort))
  table(paste(edges[, 1], edges[, 2], sep = "-"))
}

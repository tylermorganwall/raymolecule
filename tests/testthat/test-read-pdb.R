test_that("read_pdb parses a valid single-chain PDB with TER and HETATM", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(residue_1$next_serial, "A", 2, "GLY", c(1.5, 0.1, 0.2))
  residue_3 = backbone_residue_lines(residue_2$next_serial, "A", 3, "SER", c(3.0, 0.2, 0.5))

  lines = c(
    residue_1$lines,
    residue_2$lines,
    residue_3$lines,
    format_atom_line(100, "ZN", "ZN", "Z", 1, 5, 5, 5, "ZN", record = "HETATM"),
    format_conect_line(1, c(2, 3)),
    format_ter_line(101, "SER", "A", 3),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(nrow(model$atoms), 13)
  expect_equal(sum(model$atoms$record == "HETATM"), 1)
  expect_equal(nrow(model$bonds), 2)
  expect_equal(nrow(model$residues), 3)
  expect_true(all(model$residues$has_ca))
  expect_equal(model$residues$ss_class, rep("loop", 3))
  expect_true(model$residues$chain_break_before[1])
  expect_true(model$residues$chain_break_after[3])
  expect_equal(nrow(model$chains), 1)
  expect_equal(model$chains$chain_id, "A")
  expect_true(model$chains$terminated_by_ter)
  expect_identical(model$pdb_type, "pdb")
})

test_that("read_pdb parses multiple chains separated by TER", {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(chain_a_1$next_serial, "A", 2, "GLY", c(1.5, 0.1, 0.2))
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(chain_b_1$next_serial, "B", 2, "LEU", c(1.5, 4.2, 0.3))

  lines = c(
    chain_a_1$lines,
    chain_a_2$lines,
    format_ter_line(19, "GLY", "A", 2),
    chain_b_1$lines,
    chain_b_2$lines,
    format_ter_line(40, "LEU", "B", 2),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(nrow(model$chains), 2)
  expect_equal(model$chains$chain_id, c("A", "B"))
  expect_true(all(model$chains$terminated_by_ter))
  expect_true(model$residues$chain_break_after[2])
  expect_true(model$residues$chain_break_before[3])
})

test_that("read_pdb parses HELIX, SHEET, and SSBOND records", {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(chain_a_1$next_serial, "A", 2, "GLY", c(1.5, 0.1, 0.2))
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(chain_b_1$next_serial, "B", 2, "LEU", c(1.5, 4.2, 0.3))

  lines = c(
    format_helix_line(1, "H1", "ALA", "A", 1, "GLY", "A", 2, length = 2),
    format_sheet_line(1, "S1", 1, "VAL", "B", 1, "LEU", "B", 2),
    format_ssbond_line(1, "A", 1, "B", 2),
    chain_a_1$lines,
    chain_a_2$lines,
    format_ter_line(19, "GLY", "A", 2),
    chain_b_1$lines,
    chain_b_2$lines,
    format_ter_line(40, "LEU", "B", 2),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(nrow(model$helices), 1)
  expect_equal(model$helices$helix_id, "H1")
  expect_equal(nrow(model$sheets), 1)
  expect_equal(model$sheets$sheet_id, "S1")
  expect_equal(nrow(model$ssbonds), 1)
  expect_equal(model$ssbonds$chain_id_1, "A")
  expect_equal(model$ssbonds$chain_id_2, "B")
  expect_equal(model$residues$ss_class[1:2], c("helix", "helix"))
  expect_equal(model$residues$ss_class[3:4], c("sheet", "sheet"))
  expect_equal(model$residues$helix_id[1:2], c("H1", "H1"))
  expect_equal(model$residues$sheet_id[3:4], c("S1", "S1"))
})

test_that("read_pdb errors on malformed fixed-width atom lines", {
  file = write_pdb_fixture(c("ATOM      1  N  ", "END"))
  on.exit(unlink(file), add = TRUE)

  expect_error(read_pdb(file), "Malformed ATOM record")
})

test_that("read_pdb errors on multiple MODEL records", {
  residue = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  lines = c(
    format_model_line(1),
    residue$lines,
    "ENDMDL",
    format_model_line(2),
    residue$lines,
    "ENDMDL",
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  expect_error(
    read_pdb(file),
    "Multiple MODEL records are not supported yet"
  )
})

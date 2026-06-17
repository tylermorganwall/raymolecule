test_that("read_pdb parses a valid single-chain PDB with TER and HETATM", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(
    residue_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  residue_3 = backbone_residue_lines(
    residue_2$next_serial,
    "A",
    3,
    "SER",
    c(3.0, 0.2, 0.5)
  )

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
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(
    chain_b_1$next_serial,
    "B",
    2,
    "LEU",
    c(1.5, 4.2, 0.3)
  )

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
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(
    chain_b_1$next_serial,
    "B",
    2,
    "LEU",
    c(1.5, 4.2, 0.3)
  )

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

test_that("read_pdb accepts overlapping HELIX records on the same residues", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(
    residue_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  residue_3 = backbone_residue_lines(
    residue_2$next_serial,
    "A",
    3,
    "SER",
    c(3.0, 0.2, 0.5)
  )
  residue_4 = backbone_residue_lines(
    residue_3$next_serial,
    "A",
    4,
    "THR",
    c(4.5, 0.4, 0.8)
  )

  lines = c(
    format_helix_line(1, "H1", "ALA", "A", 1, "SER", "A", 3, length = 3),
    format_helix_line(2, "H2", "SER", "A", 3, "THR", "A", 4, length = 2),
    residue_1$lines,
    residue_2$lines,
    residue_3$lines,
    residue_4$lines,
    format_ter_line(50, "THR", "A", 4),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(model$residues$ss_class, rep("helix", 4))
  expect_equal(model$residues$helix_id, c("H1", "H1", "H1", "H2"))
})

test_that("read_pdb keeps the first secondary-structure assignment on overlaps", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(
    residue_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  residue_3 = backbone_residue_lines(
    residue_2$next_serial,
    "A",
    3,
    "SER",
    c(3.0, 0.2, 0.5)
  )

  lines = c(
    format_helix_line(1, "H1", "ALA", "A", 1, "GLY", "A", 2, length = 2),
    format_sheet_line(1, "S1", 1, "GLY", "A", 2, "SER", "A", 3),
    residue_1$lines,
    residue_2$lines,
    residue_3$lines,
    format_ter_line(40, "SER", "A", 3),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(model$residues$ss_class, c("helix", "helix", "sheet"))
  expect_equal(model$residues$helix_id, c("H1", "H1", NA))
  expect_equal(model$residues$sheet_id, c(NA, NA, "S1"))
})

test_that("read_pdb errors on malformed fixed-width atom lines", {
  file = write_pdb_fixture(c("ATOM      1  N  ", "END"))
  on.exit(unlink(file), add = TRUE)

  expect_error(read_pdb(file), "Malformed ATOM record")
})

test_that("read_pdb parses all MODEL records from coordinate ensembles", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(1, "A", 1, "ALA", c(10, 0, 0))
  lines = c(
    format_model_line(1),
    residue_1$lines,
    format_ter_line(5, "ALA", "A", 1),
    "ENDMDL",
    format_model_line(2),
    residue_2$lines,
    format_ter_line(5, "ALA", "A", 1),
    "ENDMDL",
    format_conect_line(1, c(2)),
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)

  expect_equal(nrow(model$atoms), 8)
  expect_equal(nrow(model$residues), 2)
  expect_equal(model$residues$model_index, c(1L, 2L))
  expect_equal(model$residues$ca_x, c(0, 10))
  expect_equal(model$atoms$serial, rep(1:4, 2))
  expect_equal(model$atoms$index, seq_len(8))
  expect_equal(nrow(model$bonds), 2)
  expect_equal(model$bonds$model_index, c(1L, 2L))
  expect_equal(model$bonds$from, c(1L, 5L))
  expect_equal(model$bonds$to, c(2L, 6L))
})

test_that("read_pdb stores metadata and reports verbose summaries", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(1, "A", 1, "ALA", c(10, 0, 0))
  header = make_pdb_line()
  header = set_pdb_field(header, 1, 6, "HEADER")
  header = set_pdb_field(header, 11, 50, "TEST STRUCTURE")
  header = set_pdb_field(header, 51, 59, "01-JAN-00")
  header = set_pdb_field(header, 63, 66, "9TST")
  lines = c(
    header,
    "TITLE     TEST PROTEIN ENSEMBLE",
    "EXPDTA    SOLUTION NMR",
    "NUMMDL    2",
    format_model_line(1),
    residue_1$lines,
    "ENDMDL",
    format_model_line(2),
    residue_2$lines,
    "ENDMDL",
    "END"
  )

  file = write_pdb_fixture(lines)
  on.exit(unlink(file), add = TRUE)

  expect_message(
    model <- read_pdb(file, verbose = TRUE),
    "Read TEST PROTEIN ENSEMBLE PDB models \\[1-2\\]"
  )
  expect_equal(model$metadata$pdb_id, "9TST")
  expect_equal(model$metadata$name, "TEST PROTEIN ENSEMBLE")
  expect_equal(model$metadata$experiment, "SOLUTION NMR")
  expect_equal(model$metadata$declared_model_count, 2L)
})

test_that("read_pdb can expand biological assemblies from REMARK 350 BIOMT", {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(
    chain_b_1$next_serial,
    "B",
    2,
    "LEU",
    c(1.5, 4.2, 0.3)
  )

  lines = c(
    format_remark350_biomolecule_line(1),
    format_remark350_apply_chains_line(c("A", "B")),
    format_remark350_biomt_line(1, 1, c(1, 0, 0, 0)),
    format_remark350_biomt_line(2, 1, c(0, 1, 0, 0)),
    format_remark350_biomt_line(3, 1, c(0, 0, 1, 0)),
    format_remark350_biomt_line(1, 2, c(1, 0, 0, 10)),
    format_remark350_biomt_line(2, 2, c(0, 1, 0, 0)),
    format_remark350_biomt_line(3, 2, c(0, 0, 1, 0)),
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

  asym_model = read_pdb(file)
  bio_model = read_pdb(file, assembly = "biological")

  expect_equal(asym_model$chains$chain_id, c("A", "B"))
  expect_equal(bio_model$chains$chain_id, c("A", "B", "A_2", "B_2"))
  expect_equal(nrow(bio_model$residues), 8)
  expect_equal(nrow(bio_model$atoms), 32)
  expect_equal(
    sort(unique(bio_model$residues$chain_id)),
    c("A", "A_2", "B", "B_2")
  )
  expect_equal(
    bio_model$residues$ca_x[bio_model$residues$chain_id == "A_2"][1],
    bio_model$residues$ca_x[bio_model$residues$chain_id == "A"][1] + 10,
    tolerance = 1e-6
  )
  expect_equal(bio_model$assembly_mode, "biological")
  expect_equal(bio_model$assembly_id, 1L)
  expect_equal(nrow(bio_model$assemblies), 4)
})

test_that("read_pdb errors when requesting a missing biological assembly", {
  residue = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  file = write_pdb_fixture(c(residue$lines, "END"))
  on.exit(unlink(file), add = TRUE)

  expect_error(
    read_pdb(file, assembly = "biological"),
    "Biological assembly transforms were not found in this PDB file"
  )
})

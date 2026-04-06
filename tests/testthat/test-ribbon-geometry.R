make_residue_table = function(ca_points, o_points, chain_id = "A") {
  count = nrow(ca_points)
  data.frame(
    chain_id = rep(chain_id, count),
    res_seq = seq_len(count),
    i_code = rep("", count),
    res_name = rep("ALA", count),
    residue_id = sprintf("%s:%d:", chain_id, seq_len(count)),
    n_x = NA_real_,
    n_y = NA_real_,
    n_z = NA_real_,
    ca_x = ca_points[, 1],
    ca_y = ca_points[, 2],
    ca_z = ca_points[, 3],
    c_x = NA_real_,
    c_y = NA_real_,
    c_z = NA_real_,
    o_x = o_points[, 1],
    o_y = o_points[, 2],
    o_z = o_points[, 3],
    has_n = FALSE,
    has_ca = TRUE,
    has_c = FALSE,
    has_o = TRUE,
    chain_break_before = c(TRUE, rep(FALSE, count - 1)),
    chain_break_after = c(rep(FALSE, count - 1), TRUE),
    ss_class = rep("loop", count),
    helix_id = rep(NA_character_, count),
    sheet_id = rep(NA_character_, count),
    stringsAsFactors = FALSE
  )
}

test_that("centripetal Catmull-Rom centerline passes through CA knot positions", {
  points = matrix(
    c(
      0, 0, 0,
      1, 0, 0,
      2, 1, 0,
      3, 1, 0
    ),
    ncol = 3,
    byrow = TRUE
  )

  chain = raymolecule:::catmull_rom_chain(points)
  sampled = raymolecule:::sample_catmull_rom_chain(chain, subdivisions = 3)

  expect_equal(sampled$centerline[c(1, 4, 7, 10), ], points, tolerance = 1e-8)
})

test_that("rotation-minimizing frame normals stay sign-consistent on a helix-like curve", {
  theta = seq(0, 2 * pi, length.out = 6)
  ca_points = cbind(cos(theta), sin(theta), seq(0, 2.5, length.out = 6))
  radial = cbind(cos(theta), sin(theta), rep(0, length(theta)))
  o_points = ca_points + 0.8 * radial
  residues = make_residue_table(ca_points, o_points)

  chain = raymolecule:::catmull_rom_chain(ca_points)
  sampled = raymolecule:::sample_catmull_rom_chain(chain, subdivisions = 4)
  guides = raymolecule:::build_residue_guides(residues, sampled$residue_parameter)
  frame = raymolecule:::build_rotation_minimizing_frame(
    centerline = sampled$centerline,
    tangent = sampled$tangent,
    guides = guides,
    chain_id = "A",
    residue_parameter = sampled$residue_parameter
  )

  dot_products = rowSums(
    frame$normal[-1, , drop = FALSE] * frame$normal[-nrow(frame$normal), , drop = FALSE]
  )
  step_angles = acos(pmin(1, pmax(-1, dot_products))) * 180 / pi

  expect_true(all(dot_products > 0))
  expect_lt(max(step_angles), 35)
})

test_that("axial guide symmetry suppresses alternating sheet-like flips", {
  ca_points = cbind(seq(0, 5, length.out = 6), rep(0, 6), rep(0, 6))
  guide_y = c(1, -1, 1, -1, 1, -1)
  o_points = ca_points + cbind(rep(0.8, 6), guide_y, rep(0, 6))
  residues = make_residue_table(ca_points, o_points)
  residues$ss_class[] = "sheet"
  residues$sheet_id[] = "AA"
  residues$sheet_strand = rep(1L, nrow(residues))

  chain = raymolecule:::catmull_rom_chain(ca_points)
  sampled = raymolecule:::sample_catmull_rom_chain(chain, subdivisions = 4)
  guides = raymolecule:::build_residue_guides(residues, sampled$residue_parameter)
  frame = raymolecule:::build_rotation_minimizing_frame(
    centerline = sampled$centerline,
    tangent = sampled$tangent,
    guides = guides,
    chain_id = "A",
    residue_parameter = sampled$residue_parameter
  )

  dot_products = rowSums(
    frame$normal[-1, , drop = FALSE] * frame$normal[-nrow(frame$normal), , drop = FALSE]
  )
  step_angles = acos(pmin(1, pmax(-1, dot_products))) * 180 / pi

  expect_true(all(dot_products > 0))
  expect_lt(max(step_angles), 25)
})

test_that("chain breaks create separate ribbon shapes", {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(chain_a_1$next_serial, "A", 2, "GLY", c(1.5, 0.1, 0.2))
  chain_b_1 = backbone_residue_lines(20, "B", 1, "VAL", c(0, 4, 0))
  chain_b_2 = backbone_residue_lines(chain_b_1$next_serial, "B", 2, "LEU", c(1.5, 4.2, 0.3))

  file = write_pdb_fixture(c(
    chain_a_1$lines,
    chain_a_2$lines,
    format_ter_line(19, "GLY", "A", 2),
    chain_b_1$lines,
    chain_b_2$lines,
    format_ter_line(40, "LEU", "B", 2),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)
  scene = generate_ribbon_scene(model, pathtrace = FALSE, subdivisions = 2)

  expect_equal(length(scene$shapes), 2)
})

test_that("ribbon sampling is refined beyond the requested minimum step count", {
  ca_points = matrix(
    c(
      0, 0, 0,
      3.8, 0, 0
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points + matrix(c(0, 1, 0, 0, 1, 0), ncol = 3, byrow = TRUE)
  residues = make_residue_table(ca_points, o_points)

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 2
  )

  expect_gt(nrow(mesh$chains[[1]]$sampled$centerline), 3)
})

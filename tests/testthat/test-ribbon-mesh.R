make_mesh_residue_table = function(ca_points, o_points, chain_id = "A") {
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

test_that("ribbon mesh is indexed, capped, watertight, and UV-bounded", {
  residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  residue_2 = backbone_residue_lines(residue_1$next_serial, "A", 2, "GLY", c(1.5, 0.1, 0.2))
  residue_3 = backbone_residue_lines(residue_2$next_serial, "A", 3, "SER", c(3.0, 0.2, 0.5))
  residue_4 = backbone_residue_lines(residue_3$next_serial, "A", 4, "THR", c(4.5, 0.4, 0.8))

  file = write_pdb_fixture(c(
    residue_1$lines,
    residue_2$lines,
    residue_3$lines,
    residue_4$lines,
    format_ter_line(30, "THR", "A", 4),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  model = read_pdb(file)
  mesh = raymolecule:::build_ribbon_mesh(
    residues = model$residues,
    subdivisions = 3
  )$chains[[1]]

  expect_false(anyNA(mesh$vertices))
  expect_false(anyNA(mesh$indices))
  expect_false(anyNA(mesh$normals))
  expect_false(anyNA(mesh$texcoords))
  expect_true(all(mesh$indices >= 0L))
  expect_true(all(mesh$indices < nrow(mesh$vertices)))

  ring_count = nrow(mesh$sampled$centerline)
  ring_size = as.integer((nrow(mesh$vertices) - 2L) / ring_count)
  expect_equal(nrow(mesh$vertices), ring_count * ring_size + 2)
  expect_equal(nrow(mesh$indices), (ring_count - 1) * ring_size * 2L + ring_size * 2L)

  edge_counts = mesh_edge_table(mesh$indices)
  expect_true(all(as.integer(edge_counts) == 2L))

  ring_u = mesh$texcoords[seq(1, ring_count * (ring_size + 1L), by = ring_size + 1L), 1]
  expect_true(all(diff(ring_u) >= -1e-8))
  expect_true(all(mesh$texcoords >= 0 - 1e-8))
  expect_true(all(mesh$texcoords <= 1 + 1e-8))
})

test_that("cross-section resolution controls ribbon ring density", {
  ca_points = matrix(
    c(
      0, 0, 0,
      1.5, 0.1, 0.2,
      3.0, 0.2, 0.4
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points + matrix(
    c(
      0, 0.8, 0,
      0, 0.8, 0,
      0, 0.8, 0
    ),
    ncol = 3,
    byrow = TRUE
  )
  residues = make_mesh_residue_table(ca_points, o_points)

  low_res_mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 2,
    cross_section_resolution = 12
  )$chains[[1]]
  high_res_mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 2,
    cross_section_resolution = 28
  )$chains[[1]]

  low_ring_count = nrow(low_res_mesh$sampled$centerline)
  high_ring_count = nrow(high_res_mesh$sampled$centerline)
  low_ring_size = as.integer((nrow(low_res_mesh$vertices) - 2L) / low_ring_count)
  high_ring_size = as.integer((nrow(high_res_mesh$vertices) - 2L) / high_ring_count)

  expect_equal(low_ring_size, 12)
  expect_equal(high_ring_size, 28)
  expect_gt(nrow(high_res_mesh$vertices), nrow(low_res_mesh$vertices))
})

test_that("cross-section resolution must be at least eight vertices", {
  ca_points = matrix(
    c(
      0, 0, 0,
      1.5, 0.1, 0.2
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points + matrix(c(0, 0.8, 0, 0, 0.8, 0), ncol = 3, byrow = TRUE)
  residues = make_mesh_residue_table(ca_points, o_points)

  expect_error(
    raymolecule:::build_ribbon_mesh(
      residues = residues,
      subdivisions = 2,
      cross_section_resolution = 6
    ),
    "cross_section_resolution must be an integer greater than or equal to 8"
  )
})

test_that("swept ribbon surface normals point outward on curved ribbons", {
  theta = seq(0, pi / 2, length.out = 5)
  ca_points = cbind(3 * cos(theta), 3 * sin(theta), seq(0, 1, length.out = 5))
  o_points = ca_points + cbind(cos(theta), sin(theta), rep(0, length(theta)))
  residues = make_mesh_residue_table(ca_points, o_points)

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 3,
    cross_section_resolution = 16
  )$chains[[1]]

  ring_count = nrow(mesh$sampled$centerline)
  ring_size = as.integer((nrow(mesh$vertices) - 2L) / ring_count)
  side_vertex_count = ring_count * ring_size
  side_vertices = mesh$vertices[seq_len(side_vertex_count), , drop = FALSE]
  side_normals = mesh$normals[seq_len(side_vertex_count), , drop = FALSE]
  centerline = mesh$sampled$centerline[rep(seq_len(ring_count), each = ring_size), , drop = FALSE]
  reference = side_vertices - centerline

  expect_true(all(rowSums(side_normals * reference) > 0))
})

test_that("warped ribbon quads use the shorter diagonal", {
  vertices = rbind(
    c(0, 0, 0),
    c(1, 0, 0),
    c(1, 1, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 1),
    c(1, 2, 1),
    c(0, 1, 1)
  )

  connectivity = raymolecule:::connect_ribbon_rings(
    vertices = vertices,
    ring_count = 2,
    ring_size = 4
  )

  expect_equal(connectivity$indices[1, ], c(0L, 1L, 4L))
  expect_equal(connectivity$indices[2, ], c(1L, 5L, 4L))
})

test_that("twist unwrapping follows full turns continuously", {
  wrapped_angles = c(170, -170, 175, -175) * pi / 180
  unwrapped = numeric(length(wrapped_angles))

  for (i in seq_along(wrapped_angles)) {
    previous_angle = if (i == 1L) NULL else unwrapped[i - 1L]
    unwrapped[i] = raymolecule:::choose_continuous_twist_angle(
      angle = wrapped_angles[i],
      previous_angle = previous_angle
    )
  }

  expect_lt(max(abs(diff(unwrapped))) * 180 / pi, 30)
})

test_that("twist interpolation is smooth at interior knots", {
  positions = seq(0, 4, by = 0.01)
  knot_positions = 0:4
  knot_angles = c(0, 0.6, -0.2, 0.5, 0.1)

  interpolated = raymolecule:::interpolate_twist_angles(
    positions = positions,
    knot_positions = knot_positions,
    knot_angles = knot_angles
  )

  left_index = which.min(abs(positions - 1.99))
  knot_index = which.min(abs(positions - 2.00))
  right_index = which.min(abs(positions - 2.01))

  left_slope = (interpolated[knot_index] - interpolated[left_index]) / 0.01
  right_slope = (interpolated[right_index] - interpolated[knot_index]) / 0.01

  expect_lt(abs(left_slope - right_slope), 0.2)
})

test_that("global twist smoothing reduces local roll jumps", {
  positions = seq(0, 11)
  raw_angles = seq(0, pi / 2, length.out = 12) +
    c(0, 0.7, -0.6, 0.5, -0.4, 0.3, -0.3, 0.4, -0.5, 0.6, -0.7, 0)

  smoothed_angles = raymolecule:::smooth_twist_angles(
    positions = positions,
    angles = raw_angles
  )

  expect_lt(max(abs(diff(smoothed_angles))), max(abs(diff(raw_angles))))
  expect_equal(smoothed_angles[c(1, length(smoothed_angles))], raw_angles[c(1, length(raw_angles))], tolerance = 0.25)
})

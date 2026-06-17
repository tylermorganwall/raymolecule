make_mesh_residue_table = function(
  ca_points,
  o_points,
  chain_id = "A",
  ss_class = NULL,
  sheet_id = NULL,
  sheet_strand = NULL
) {
  count = nrow(ca_points)
  if (is.null(ss_class)) {
    ss_class = rep("loop", count)
  }
  if (is.null(sheet_id)) {
    sheet_id = rep(NA_character_, count)
  }
  if (is.null(sheet_strand)) {
    sheet_strand = rep(NA_integer_, count)
  }
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
    ss_class = ss_class,
    helix_id = rep(NA_character_, count),
    sheet_id = sheet_id,
    sheet_strand = sheet_strand,
    stringsAsFactors = FALSE
  )
}

ribbon_step_metrics = function(mesh) {
  sampled = mesh$sampled
  centerline = sampled$centerline
  step_vectors = centerline[-1, , drop = FALSE] -
    centerline[-nrow(centerline), , drop = FALSE]
  step_lengths = sqrt(rowSums(step_vectors^2))
  tangent = sampled$tangent
  alignment = rowSums(
    step_vectors * tangent[-nrow(tangent), , drop = FALSE]
  ) /
    pmax(step_lengths, 1e-12)
  tangent_turn = rowSums(
    tangent[-1, , drop = FALSE] * tangent[-nrow(tangent), , drop = FALSE]
  )

  list(
    length = step_lengths,
    alignment = alignment,
    tangent_turn = tangent_turn,
    profile_kind = sampled$sampled$profile_kind
  )
}

build_transition_mesh = function(ca_points, ss_class) {
  o_points = ca_points +
    matrix(
      rep(c(0.4, 0.8, 0.2), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points = ca_points,
    o_points = o_points,
    ss_class = ss_class
  )

  raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 16,
    cross_section_resolution = 16
  )$chains[[1]]
}

test_that("ribbon mesh is indexed, capped, watertight, and UV-bounded", {
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
  expect_equal(
    nrow(mesh$indices),
    (ring_count - 1) * ring_size * 2L + ring_size * 2L
  )

  edge_counts = mesh_edge_table(mesh$indices)
  expect_true(all(as.integer(edge_counts) == 2L))

  ring_u = mesh$texcoords[
    seq(1, ring_count * (ring_size + 1L), by = ring_size + 1L),
    1
  ]
  expect_true(all(diff(ring_u) >= -1e-8))
  expect_true(all(mesh$texcoords >= 0 - 1e-8))
  expect_true(all(mesh$texcoords <= 1 + 1e-8))
})

test_that("cross-section resolution controls ribbon ring density", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0.1,
      0.2,
      3.0,
      0.2,
      0.4
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      c(
        0,
        0.8,
        0,
        0,
        0.8,
        0,
        0,
        0.8,
        0
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
  low_ring_size = as.integer(
    (nrow(low_res_mesh$vertices) - 2L) / low_ring_count
  )
  high_ring_size = as.integer(
    (nrow(high_res_mesh$vertices) - 2L) / high_ring_count
  )

  expect_equal(low_ring_size, 12)
  expect_equal(high_ring_size, 28)
  expect_gt(nrow(high_res_mesh$vertices), nrow(low_res_mesh$vertices))
})

test_that("cross-section resolution must be at least eight vertices", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0.1,
      0.2
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

test_that("profile generators share a common perimeter phase", {
  ellipse = raymolecule:::ellipse_profile_2d(24, 0.6, 0.6)
  rectangle = raymolecule:::rectangle_profile_2d(24, 1.6, 0.44)
  rounded = raymolecule:::rounded_rectangle_profile_2d(24, 1.6, 0.25)

  expect_equal(which.max(ellipse[, 1] + ellipse[, 2]), 1L)
  expect_equal(which.max(rectangle[, 1] + rectangle[, 2]), 1L)
  expect_equal(which.max(rounded[, 1] + rounded[, 2]), 1L)
})

test_that("helix centerlines preserve the CA helix radius", {
  radius = 2.3
  theta = seq(0, by = 100 * pi / 180, length.out = 14)
  ca_points = cbind(
    radius * cos(theta),
    radius * sin(theta),
    seq(0, by = 1.5, length.out = length(theta))
  )
  radial_guides = cbind(cos(theta), sin(theta), rep(0, length(theta)))
  o_points = ca_points + 0.8 * radial_guides
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = rep("helix", nrow(ca_points))
  )

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 8,
    cross_section_resolution = 16
  )$chains[[1]]

  center_radius = sqrt(rowSums(mesh$sampled$centerline[, 1:2]^2))
  midpoint_points = (ca_points[-nrow(ca_points), ] + ca_points[-1, ]) / 2
  midpoint_radius = sqrt(rowSums(midpoint_points[, 1:2]^2))

  expect_gt(min(center_radius), radius * 0.98)
  expect_lt(stats::sd(center_radius), radius * 0.02)
  expect_gt(mean(center_radius), radius * 0.85)
  expect_gt(mean(center_radius), mean(midpoint_radius) * 1.25)
})

test_that("sheet-to-helix transitions keep moving forward", {
  ca_points = matrix(
    c(
      -0.320,
      36.894,
      5.259,
      1.849,
      33.903,
      4.222,
      5.202,
      33.307,
      5.993,
      7.684,
      30.754,
      4.629,
      8.531,
      28.869,
      7.801,
      7.400,
      28.022,
      11.358,
      9.523,
      27.946,
      14.533,
      9.128,
      24.164,
      14.895,
      10.772,
      23.460,
      11.551,
      13.392,
      26.130,
      12.283,
      14.102,
      24.625,
      15.735
    ),
    ncol = 3,
    byrow = TRUE
  )
  mesh = build_transition_mesh(
    ca_points = ca_points,
    ss_class = c(
      "sheet",
      "loop",
      "loop",
      "loop",
      "sheet",
      "sheet",
      rep("helix", 5)
    )
  )
  metrics = ribbon_step_metrics(mesh)
  transition_rows = which(
    metrics$profile_kind[-length(metrics$profile_kind)] == "strand" &
      metrics$profile_kind[-1] == "helix"
  )

  expect_length(transition_rows, 1)
  expect_gt(metrics$alignment[transition_rows], 0.95)
  expect_gt(metrics$tangent_turn[transition_rows], 0.95)
  expect_lt(metrics$length[transition_rows], 0.4)
})

test_that("helix-to-coil transitions keep tangents aligned with the blended path", {
  ca_points = matrix(
    c(
      8.129,
      30.791,
      28.778,
      7.929,
      27.119,
      27.774,
      8.819,
      28.000,
      24.158,
      11.824,
      30.074,
      25.357,
      13.123,
      27.167,
      27.470,
      12.747,
      24.836,
      24.462,
      14.701,
      27.168,
      22.114,
      17.502,
      27.340,
      24.665,
      17.720,
      23.644,
      25.605,
      16.532,
      22.084,
      22.308,
      14.140,
      19.841,
      24.270,
      10.499,
      20.550,
      25.159,
      10.217,
      21.653,
      28.777
    ),
    ncol = 3,
    byrow = TRUE
  )
  mesh = build_transition_mesh(
    ca_points = ca_points,
    ss_class = c(rep("helix", 10), rep("loop", 3))
  )
  metrics = ribbon_step_metrics(mesh)
  transition_rows = which(
    metrics$profile_kind[-length(metrics$profile_kind)] == "helix" &
      metrics$profile_kind[-1] == "coil"
  )

  expect_length(transition_rows, 1)
  expect_gt(metrics$alignment[transition_rows], 0.95)
  expect_gt(metrics$tangent_turn[transition_rows], 0.95)
  expect_lt(metrics$length[transition_rows], 0.4)
})

test_that("sheet runs generate tapered arrow widths", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3,
      6.0,
      0.0,
      0.4,
      7.5,
      -0.1,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = c("loop", "sheet", "sheet", "sheet", "loop", "loop"),
    sheet_id = c(NA, "S1", "S1", "S1", NA, NA)
  )

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 4,
    cross_section_resolution = 16
  )$chains[[1]]

  widths = mesh$sampled$sampled$ribbon_width
  parameters = mesh$sampled$residue_parameter
  head_start_index = which(parameters > 2.5)[1]
  head_mask = parameters >= parameters[head_start_index] & parameters <= 3.25
  taper_mask = parameters >= 2.75 & parameters <= 3.25

  expect_gt(max(widths[head_mask]), 1.6)
  expect_equal(widths[which.min(abs(parameters - 2.5))], 1.6, tolerance = 1e-6)
  expect_lt(parameters[head_start_index] - 2.5, 0.01)
  expect_gt(widths[head_start_index], 1.6)
  expect_true(all(diff(widths[taper_mask]) < 0))
  expect_lt(widths[which.min(abs(parameters - 3.25))], 0.7)
})

test_that("peptide-plane transitions avoid collapsed interior rings", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3,
      6.0,
      0.0,
      0.4,
      7.5,
      -0.1,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = c("loop", "sheet", "sheet", "sheet", "loop", "loop"),
    sheet_id = c(NA, "S1", "S1", "S1", NA, NA)
  )

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 8,
    cross_section_resolution = 24
  )$chains[[1]]

  expect_true(all(mesh$sampled$sampled$ribbon_width > 0))
  expect_true(all(mesh$sampled$sampled$ribbon_thickness > 0))

  index_rows = mesh$indices + 1L
  edge1 = mesh$vertices[index_rows[, 2], , drop = FALSE] -
    mesh$vertices[index_rows[, 1], , drop = FALSE]
  edge2 = mesh$vertices[index_rows[, 3], , drop = FALSE] -
    mesh$vertices[index_rows[, 1], , drop = FALSE]
  face_normals = cbind(
    edge1[, 2] * edge2[, 3] - edge1[, 3] * edge2[, 2],
    edge1[, 3] * edge2[, 1] - edge1[, 1] * edge2[, 3],
    edge1[, 1] * edge2[, 2] - edge1[, 2] * edge2[, 1]
  )
  face_area = sqrt(rowSums(face_normals^2))

  expect_true(all(face_area > 1e-8))
})

test_that("peptide-plane transition side faces keep outward winding", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3,
      6.0,
      0.0,
      0.4,
      7.5,
      -0.1,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = c("loop", "sheet", "sheet", "sheet", "loop", "loop"),
    sheet_id = c(NA, "S1", "S1", "S1", NA, NA),
    sheet_strand = c(NA, 1L, 1L, 1L, NA, NA)
  )

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 8,
    cross_section_resolution = 60
  )$chains[[1]]

  ring_count = nrow(mesh$sampled$centerline)
  ring_size = as.integer((nrow(mesh$vertices) - 2L) / ring_count)
  side_triangle_count = (ring_count - 1L) * ring_size * 2L
  side_indices = mesh$indices[seq_len(side_triangle_count), , drop = FALSE] + 1L
  edge1 = mesh$vertices[side_indices[, 2], , drop = FALSE] -
    mesh$vertices[side_indices[, 1], , drop = FALSE]
  edge2 = mesh$vertices[side_indices[, 3], , drop = FALSE] -
    mesh$vertices[side_indices[, 1], , drop = FALSE]
  face_normals = cbind(
    edge1[, 2] * edge2[, 3] - edge1[, 3] * edge2[, 2],
    edge1[, 3] * edge2[, 1] - edge1[, 1] * edge2[, 3],
    edge1[, 1] * edge2[, 2] - edge1[, 2] * edge2[, 1]
  )
  reference_normals =
    mesh$normals[side_indices[, 1], , drop = FALSE] +
    mesh$normals[side_indices[, 2], , drop = FALSE] +
    mesh$normals[side_indices[, 3], , drop = FALSE]

  expect_true(all(rowSums(face_normals * reference_normals) >= -1e-8))
})

test_that("strand-exit segment boundaries use matching profile families", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3,
      6.0,
      0.0,
      0.4,
      7.5,
      -0.1,
      0.5
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = c("loop", "sheet", "sheet", "sheet", "loop", "loop"),
    sheet_id = c(NA, "S1", "S1", "S1", NA, NA)
  )

  planes = raymolecule:::build_peptide_planes(residues)
  exit_profiles = raymolecule:::ribbon_segment_profiles(
    planes = planes,
    plane_index = 3L,
    cross_section_resolution = 24,
    ribbon_width = 1.6,
    ribbon_thickness = 0.25
  )
  next_profiles = raymolecule:::ribbon_segment_profiles(
    planes = planes,
    plane_index = 4L,
    cross_section_resolution = 24,
    ribbon_width = 1.6,
    ribbon_thickness = 0.25
  )

  expect_equal(next_profiles$start, exit_profiles$tip, tolerance = 1e-8)
})

test_that("complete backbone segments default to the peptide-plane ribbon builder", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(ca_points, o_points)

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 4,
    cross_section_resolution = 16
  )$chains[[1]]

  expect_true("profile_kind" %in% colnames(mesh$sampled$sampled))
  expect_true(all(mesh$sampled$sampled$profile_kind == "coil"))
})

test_that("CA/RMF ribbons remain the fallback for incomplete peptide geometry", {
  ca_points = matrix(
    c(
      0,
      0,
      0,
      1.5,
      0,
      0.1,
      3.0,
      0.1,
      0.2,
      4.5,
      0.1,
      0.3
    ),
    ncol = 3,
    byrow = TRUE
  )
  o_points = ca_points +
    matrix(
      rep(c(0, 0.8, 0), nrow(ca_points)),
      ncol = 3,
      byrow = TRUE
    )
  residues = make_mesh_residue_table(
    ca_points,
    o_points,
    ss_class = c("loop", "sheet", "sheet", "sheet"),
    sheet_id = c(NA, "S1", "S1", "S1"),
    sheet_strand = c(NA, 1L, 1L, 1L)
  )
  residues$has_o[2] = FALSE
  residues$o_x[2] = NA_real_
  residues$o_y[2] = NA_real_
  residues$o_z[2] = NA_real_

  mesh = raymolecule:::build_ribbon_mesh(
    residues = residues,
    subdivisions = 4,
    cross_section_resolution = 16
  )$chains[[1]]

  expect_true(all(mesh$sampled$sampled$ribbon_width == 1.6))
  expect_true(all(mesh$sampled$sampled$ribbon_thickness == 0.25))
  expect_true(all(mesh$sampled$sampled$ribbon_profile_exponent == 4))
  expect_false("profile_kind" %in% colnames(mesh$sampled$sampled))
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
  centerline = mesh$sampled$centerline[
    rep(seq_len(ring_count), each = ring_size),
    ,
    drop = FALSE
  ]
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

test_that("connect_ribbon_rings uses normals to keep local quad winding outward", {
  vertices = structure(
    c(
      -0.187936,
      0.055093,
      -0.250689,
      0.478584,
      -0.231355,
      0.010211,
      -0.312415,
      0.516418,
      0.098852,
      -0.246141,
      0.146229,
      0.221497,
      -0.034936,
      -0.180788,
      -0.039402,
      0.187857,
      0,
      0,
      0,
      0,
      1.056609,
      1.02,
      1.120628,
      0.991434
    ),
    dim = c(8L, 3L)
  )
  normals = structure(
    c(
      -0.979313,
      0.103442,
      -0.949061,
      0.93913,
      -0.995977,
      0.086889,
      -0.970275,
      0.926099,
      0.202353,
      -0.994636,
      0.315093,
      0.343563,
      -0.079472,
      -0.982819,
      -0.071118,
      0.364025,
      0,
      0,
      0,
      0,
      0.041412,
      -0.162839,
      0.231321,
      -0.099126
    ),
    dim = c(8L, 3L)
  )

  count_bad_faces = function(connectivity) {
    index_rows = connectivity$indices + 1L
    edge1 = vertices[index_rows[, 2], , drop = FALSE] -
      vertices[index_rows[, 1], , drop = FALSE]
    edge2 = vertices[index_rows[, 3], , drop = FALSE] -
      vertices[index_rows[, 1], , drop = FALSE]
    face_normals = cbind(
      edge1[, 2] * edge2[, 3] - edge1[, 3] * edge2[, 2],
      edge1[, 3] * edge2[, 1] - edge1[, 1] * edge2[, 3],
      edge1[, 1] * edge2[, 2] - edge1[, 2] * edge2[, 1]
    )
    reference_normals =
      normals[index_rows[, 1], , drop = FALSE] +
      normals[index_rows[, 2], , drop = FALSE] +
      normals[index_rows[, 3], , drop = FALSE]
    sum(rowSums(face_normals * reference_normals) < 0)
  }

  legacy = raymolecule:::connect_ribbon_rings(
    vertices = vertices,
    ring_count = 2L,
    ring_size = 4L
  )
  guided = raymolecule:::connect_ribbon_rings(
    vertices = vertices,
    ring_count = 2L,
    ring_size = 4L,
    normals = normals
  )

  expect_gt(count_bad_faces(legacy), 0)
  expect_equal(count_bad_faces(guided), 0)
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

test_that("axial twist unwrapping suppresses guide sign flips", {
  wrapped_angles = c(5, 175, 10, -170) * pi / 180
  unwrapped = numeric(length(wrapped_angles))

  for (i in seq_along(wrapped_angles)) {
    previous_angle = if (i == 1L) NULL else unwrapped[i - 1L]
    unwrapped[i] = raymolecule:::choose_continuous_axial_twist_angle(
      angle = wrapped_angles[i],
      previous_angle = previous_angle
    )
  }

  expect_lt(max(abs(diff(unwrapped))) * 180 / pi, 30)
})

test_that("sheet twist knot smoothing reduces strand-scale wiggle", {
  residues = make_mesh_residue_table(
    ca_points = cbind(seq_len(6), rep(0, 6), rep(0, 6)),
    o_points = cbind(seq_len(6), rep(1, 6), rep(0, 6)),
    ss_class = rep("sheet", 6),
    sheet_id = rep("S1", 6),
    sheet_strand = rep(1L, 6)
  )
  knot_positions = seq(0, 5)
  knot_angles = c(0, 0.9, 0.1, 1.0, 0.2, 1.1)

  smoothed = raymolecule:::smooth_sheet_twist_knots(
    knot_positions = knot_positions,
    knot_angles = knot_angles,
    knot_residue_index = seq_len(6),
    residues = residues
  )

  expect_lt(max(abs(diff(smoothed))), max(abs(diff(knot_angles))))
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
  expect_equal(
    smoothed_angles[c(1, length(smoothed_angles))],
    raw_angles[c(1, length(raw_angles))],
    tolerance = 0.25
  )
})

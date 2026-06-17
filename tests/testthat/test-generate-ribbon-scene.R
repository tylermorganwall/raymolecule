make_two_chain_model = function() {
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

  return(read_pdb(file))
}

make_single_chain_model = function() {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )

  file = write_pdb_fixture(c(
    chain_a_1$lines,
    chain_a_2$lines,
    format_ter_line(19, "GLY", "A", 2),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file))
}

make_two_model_ribbon_model = function() {
  model_1_residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  model_1_residue_2 = backbone_residue_lines(
    model_1_residue_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )
  model_2_residue_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  model_2_residue_2 = backbone_residue_lines(
    model_2_residue_1$next_serial,
    "A",
    2,
    "GLY",
    c(0.1, 3.0, 0.2)
  )

  file = write_pdb_fixture(c(
    "TITLE     TEST RIBBON ENSEMBLE",
    "NUMMDL    2",
    format_model_line(1),
    model_1_residue_1$lines,
    model_1_residue_2$lines,
    "ENDMDL",
    format_model_line(2),
    model_2_residue_1$lines,
    model_2_residue_2$lines,
    "ENDMDL",
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file))
}

make_texture_file = function() {
  file = tempfile(fileext = ".png")
  grDevices::png(file, width = 16, height = 4)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::rect(
    xleft = c(0, 0.25, 0.5, 0.75),
    ybottom = 0,
    xright = c(0.25, 0.5, 0.75, 1),
    ytop = 1,
    col = c("red", "yellow", "green", "blue"),
    border = NA
  )
  grDevices::dev.off()
  return(file)
}

make_sheet_model = function() {
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

  residues = data.frame(
    chain_id = rep("A", nrow(ca_points)),
    res_seq = seq_len(nrow(ca_points)),
    i_code = rep("", nrow(ca_points)),
    res_name = rep("ALA", nrow(ca_points)),
    residue_id = sprintf("A:%d:", seq_len(nrow(ca_points))),
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
    chain_break_before = c(TRUE, rep(FALSE, nrow(ca_points) - 1)),
    chain_break_after = c(rep(FALSE, nrow(ca_points) - 1), TRUE),
    ss_class = c("loop", "sheet", "sheet", "sheet", "loop", "loop"),
    helix_id = rep(NA_character_, nrow(ca_points)),
    sheet_id = c(NA, "S1", "S1", "S1", NA, NA),
    sheet_strand = c(NA, 1L, 1L, 1L, NA, NA),
    stringsAsFactors = FALSE
  )

  return(list(pdb_type = "pdb", residues = residues))
}

make_biological_assembly_model = function() {
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

  file = write_pdb_fixture(c(
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
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file, assembly = "biological"))
}

make_ribbon_hetero_model = function() {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )

  file = write_pdb_fixture(c(
    chain_a_1$lines,
    chain_a_2$lines,
    format_atom_line(
      50,
      "ZN",
      "ZN",
      "A",
      101,
      0.8,
      1.2,
      0.1,
      "ZN",
      record = "HETATM"
    ),
    format_atom_line(
      51,
      "O",
      "HOH",
      "A",
      201,
      2.0,
      -1.0,
      0.0,
      "O",
      record = "HETATM"
    ),
    format_ter_line(19, "GLY", "A", 2),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file))
}

make_ribbon_phosphate_model = function() {
  chain_a_1 = backbone_residue_lines(1, "A", 1, "ALA", c(0, 0, 0))
  chain_a_2 = backbone_residue_lines(
    chain_a_1$next_serial,
    "A",
    2,
    "GLY",
    c(1.5, 0.1, 0.2)
  )

  file = write_pdb_fixture(c(
    chain_a_1$lines,
    chain_a_2$lines,
    format_atom_line(
      60,
      "P",
      "PO4",
      "A",
      201,
      2.5,
      1.8,
      0.0,
      "P",
      record = "HETATM"
    ),
    format_atom_line(
      61,
      "O1",
      "PO4",
      "A",
      201,
      3.4,
      1.8,
      0.0,
      "O",
      record = "HETATM"
    ),
    format_atom_line(
      62,
      "O2",
      "PO4",
      "A",
      201,
      2.1,
      2.6,
      0.2,
      "O",
      record = "HETATM"
    ),
    format_atom_line(
      63,
      "O3",
      "PO4",
      "A",
      201,
      2.1,
      1.0,
      -0.2,
      "O",
      record = "HETATM"
    ),
    format_atom_line(
      64,
      "O4",
      "PO4",
      "A",
      201,
      1.8,
      1.9,
      0.8,
      "O",
      record = "HETATM"
    ),
    format_conect_line(60, c(61, 62, 63, 64)),
    format_conect_line(61, 60),
    format_conect_line(62, 60),
    format_conect_line(63, 60),
    format_conect_line(64, 60),
    format_ter_line(19, "GLY", "A", 2),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file))
}

make_ribbon_protein_atom_model = function() {
  file = write_pdb_fixture(c(
    format_atom_line(1, "N", "ALA", "A", 1, 0.00, 0.00, 0.00, "N"),
    format_atom_line(2, "CA", "ALA", "A", 1, 1.45, 0.00, 0.00, "C"),
    format_atom_line(3, "C", "ALA", "A", 1, 2.10, 1.35, 0.00, "C"),
    format_atom_line(4, "O", "ALA", "A", 1, 1.60, 2.40, 0.00, "O"),
    format_atom_line(5, "N", "GLY", "A", 2, 3.40, 1.35, 0.00, "N"),
    format_atom_line(6, "CA", "GLY", "A", 2, 4.10, 2.60, 0.00, "C"),
    format_atom_line(7, "C", "GLY", "A", 2, 5.50, 2.40, 0.00, "C"),
    format_atom_line(8, "O", "GLY", "A", 2, 6.00, 3.50, 0.00, "O"),
    format_ter_line(9, "GLY", "A", 2),
    "END"
  ))
  on.exit(unlink(file), add = TRUE)

  return(read_pdb(file))
}

test_that("chain mode assigns distinct materials per chain", {
  model = make_two_chain_model()
  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2,
    color_mode = "chain"
  )

  material_a = scene$materials[[1]][[1]]$diffuse
  material_b = scene$materials[[2]][[1]]$diffuse

  expect_false(isTRUE(all.equal(material_a, material_b)))
})

test_that("biological assembly chain mode assigns materials per transformed chain", {
  model = make_biological_assembly_model()

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2,
    color_mode = "chain"
  )

  expect_equal(length(scene$materials), 4)
  diffuse_values = vapply(
    scene$materials,
    function(x) paste(round(x[[1]]$diffuse, 6), collapse = ","),
    character(1)
  )
  expect_equal(length(unique(diffuse_values)), 4)
})

test_that("ribbon scenes show non-water hetero atoms by default", {
  model = make_ribbon_hetero_model()

  default_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2
  )
  water_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    show_waters = TRUE
  )
  hidden_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    show_hetero_atoms = FALSE
  )

  expect_equal(sum(default_scene$shape == "sphere"), 1)
  expect_equal(sum(water_scene$shape == "sphere"), 2)
  expect_equal(sum(hidden_scene$shape == "sphere"), 0)
})

test_that("ribbon hetero atoms use reduced sphere radii", {
  expect_lt(
    raymolecule:::ribbon_display_atom_radius("P"),
    raymolecule:::display_atom_radius("P")
  )
  expect_lt(
    raymolecule:::ribbon_display_atom_radius("Zn"),
    raymolecule:::display_atom_radius("Zn")
  )
  expect_lt(
    raymolecule:::ribbon_display_bond_radius("P", "O"),
    raymolecule:::ribbon_display_atom_radius("P")
  )
  expect_lt(
    raymolecule:::ribbon_display_bond_radius("P", "O"),
    raymolecule:::ribbon_display_atom_radius("O")
  )
  expect_lt(
    raymolecule:::ribbon_display_atom_radius("C", scale = 0.18),
    raymolecule:::ribbon_display_atom_radius("C", scale = 0.45)
  )
  expect_lt(
    raymolecule:::ribbon_display_bond_radius("C", "N", scale = 0.14),
    raymolecule:::ribbon_display_bond_radius("C", "N", scale = 0.6)
  )
})

test_that("ribbon scenes show hetero bonds between displayed hetero atoms", {
  model = make_ribbon_phosphate_model()

  visible_atoms = raymolecule:::select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = TRUE,
    show_waters = FALSE
  )
  visible_bonds = raymolecule:::select_ribbon_display_bonds(
    model = model,
    atoms = visible_atoms,
    show_hetero_bonds = TRUE
  )
  default_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2
  )
  hidden_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    show_hetero_bonds = FALSE
  )

  expect_equal(nrow(visible_atoms), 5)
  expect_equal(nrow(visible_bonds), 4)
  expect_equal(sum(default_scene$shape == "sphere"), 5)
  expect_equal(sum(default_scene$shape == "cylinder"), 8)
  expect_equal(sum(hidden_scene$shape == "sphere"), 5)
  expect_equal(sum(hidden_scene$shape == "cylinder"), 0)
})

test_that("ribbon scenes can show protein atoms as an overlay", {
  model = make_ribbon_protein_atom_model()

  visible_atoms = raymolecule:::select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = FALSE,
    show_waters = FALSE,
    show_protein_atoms = TRUE
  )
  default_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2
  )
  atom_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    show_protein_atoms = TRUE
  )

  expect_equal(nrow(visible_atoms), 8)
  expect_true(all(visible_atoms$record == "ATOM"))
  expect_equal(sum(default_scene$shape == "sphere"), 0)
  expect_equal(sum(atom_scene$shape == "sphere"), 8)
})

test_that("ribbon scenes infer protein bond overlays from atoms", {
  model = make_ribbon_protein_atom_model()

  bond_atoms = raymolecule:::select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = FALSE,
    show_waters = FALSE,
    show_protein_atoms = TRUE
  )
  visible_bonds = raymolecule:::select_ribbon_display_bonds(
    model = model,
    atoms = bond_atoms,
    show_hetero_bonds = FALSE,
    show_protein_bonds = TRUE
  )
  default_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2
  )
  bond_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    show_protein_bonds = TRUE
  )

  expect_equal(nrow(visible_bonds), 7)
  expect_equal(sum(default_scene$shape == "cylinder"), 0)
  expect_equal(sum(bond_scene$shape == "cylinder"), 14)
  expect_equal(sum(bond_scene$shape == "sphere"), 0)
})

test_that("uv mode attaches texture coordinates and texture locations", {
  model = make_two_chain_model()
  texture = make_texture_file()
  on.exit(unlink(texture), add = TRUE)

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2,
    color_mode = "uv",
    texture = texture
  )

  expect_true(all(
    vapply(scene$materials, function(x) x[[1]]$diffuse_texname, character(1)) ==
      texture
  ))
  expect_true(length(scene$texcoords) > 0)
})

test_that("ribbon scenes omit explicit vertex normals by default", {
  model = make_two_chain_model()

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2
  )

  expect_true(all(unlist(lapply(
    scene$shapes,
    function(x) !x$has_vertex_normals
  ))))
  expect_true(all(vapply(scene$normals, nrow, integer(1)) == 0L))
})

test_that("ribbon scenes include explicit vertex normals when requested", {
  model = make_two_chain_model()

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2,
    use_vertex_normals = TRUE
  )

  expect_true(all(unlist(lapply(
    scene$shapes,
    function(x) x$has_vertex_normals
  ))))
  expect_true(all(vapply(scene$normals, nrow, integer(1)) > 0L))
})

test_that("single-chain default mode uses the built-in UV texture", {
  model = make_single_chain_model()

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2
  )

  texture_names = vapply(
    scene$materials,
    function(x) x[[1]]$diffuse_texname,
    character(1)
  )

  expect_equal(length(texture_names), 1)
  expect_true(nzchar(texture_names[[1]]))
  expect_true(file.exists(texture_names[[1]]))
})

test_that("multi-chain default mode preserves per-chain coloring", {
  model = make_two_chain_model()

  scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2
  )

  diffuse_texnames = vapply(
    scene$materials,
    function(x) x[[1]]$diffuse_texname,
    character(1)
  )
  material_a = scene$materials[[1]][[1]]$diffuse
  material_b = scene$materials[[2]][[1]]$diffuse

  expect_true(all(diffuse_texnames == ""))
  expect_false(isTRUE(all.equal(material_a, material_b)))
})

test_that("generate_ribbon_scene returns valid rayvertex and rayrender scenes", {
  model = make_two_chain_model()

  raster_scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2
  )
  path_scene = generate_ribbon_scene(model, pathtrace = TRUE, subdivisions = 2)
  smooth_path_scene = generate_ribbon_scene(
    model,
    pathtrace = TRUE,
    subdivisions = 2,
    use_vertex_normals = TRUE
  )

  expect_s3_class(raster_scene, "ray_mesh")
  expect_s3_class(path_scene, "ray_scene")
  expect_s3_class(smooth_path_scene, "ray_scene")
})

test_that("generate_ribbon_scene can render a selected PDB model", {
  model = make_two_model_ribbon_model()

  selected_model = raymolecule:::select_pdb_models(model, model_id = 2)
  expect_equal(unique(selected_model$residues$model_id), 2L)
  expect_equal(nrow(selected_model$residues), 2)

  expect_message(
    selected_scene <- generate_ribbon_scene(
      model,
      model_id = 2,
      verbose = TRUE,
      pathtrace = FALSE,
      cross_section_resolution = 8,
      subdivisions = 1
    ),
    "Rendering TEST RIBBON ENSEMBLE PDB models \\[2\\]"
  )
  expect_s3_class(selected_scene, "ray_mesh")
  expect_error(
    generate_ribbon_scene(model, model_id = 3),
    "PDB model_id value\\(s\\) not found"
  )
})

test_that("flat ribbon normals fall back cleanly on degenerate triangles", {
  vertices = matrix(
    c(
      0,
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      0
    ),
    ncol = 3,
    byrow = TRUE
  )
  indices = matrix(c(0L, 1L, 2L), ncol = 3, byrow = TRUE)
  fallback_normals = matrix(
    c(
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      1
    ),
    ncol = 3,
    byrow = TRUE
  )
  fallback_norm_indices = matrix(c(0L, 1L, 2L), ncol = 3, byrow = TRUE)

  flat_normals = raymolecule:::build_flat_mesh_normals(
    vertices = vertices,
    indices = indices,
    fallback_normals = fallback_normals,
    fallback_norm_indices = fallback_norm_indices
  )

  expect_equal(flat_normals$normals, matrix(c(0, 0, 1), ncol = 3))
  expect_equal(flat_normals$norm_indices, matrix(c(0L, 0L, 0L), ncol = 3))
})

test_that("pathtraced sheet ribbons tolerate low-resolution degenerate faces", {
  model = make_sheet_model()

  scene = expect_no_error(
    generate_ribbon_scene(
      model,
      pathtrace = TRUE,
      cross_section_resolution = 8,
      subdivisions = 4
    )
  )
  expect_s3_class(scene, "ray_scene")
})

test_that("render_model works for pathtraced ribbon scenes", {
  model = make_two_chain_model()
  path_scene = generate_ribbon_scene(model, pathtrace = TRUE, subdivisions = 2)

  path_image = render_model(
    path_scene,
    width = 32,
    height = 32,
    samples = 1,
    min_variance = 0,
    sample_method = "sobol_blue",
    parallel = FALSE
  )

  expect_true(is.array(path_image))
})

test_that("render_model forwards dimensions and explicit camera values", {
  model = make_two_chain_model()
  path_scene = generate_ribbon_scene(model, pathtrace = TRUE, subdivisions = 2)
  raster_scene = generate_ribbon_scene(
    model,
    pathtrace = FALSE,
    subdivisions = 2
  )
  pathtrace_args = NULL
  raster_args = NULL
  light_matrix = matrix(c(1, 1, 1), nrow = 1)
  lookat = c(1, 2, 3)
  lookfrom = c(4, 5, 6)

  local_mocked_bindings(
    render_scene = function(...) {
      pathtrace_args <<- list(...)
      array(0, dim = c(1, 1, 4))
    },
    rasterize_scene = function(...) {
      raster_args <<- list(...)
      array(0, dim = c(1, 1, 4))
    },
    rotate_mesh = function(scene, angle, order_rotation) {
      scene
    },
    group_objects = function(scene, angle, order_rotation) {
      scene
    },
    .package = "raymolecule"
  )

  expect_no_error(
    render_model(
      path_scene,
      lights = "none",
      lookat = lookat,
      lookfrom = lookfrom,
      samples = 1,
      min_variance = 0,
      parallel = FALSE
    )
  )
  expect_equal(pathtrace_args$lookat, lookat)
  expect_equal(pathtrace_args$lookfrom, lookfrom)
  expect_equal(pathtrace_args$width, 800L)
  expect_equal(pathtrace_args$height, 800L)
  expect_error(
    render_model(path_scene, width = 0),
    "width must be a positive finite number"
  )

  expect_no_error(
    render_model(
      raster_scene,
      width = 64,
      height = 65,
      lights = light_matrix,
      lookat = lookat,
      lookfrom = lookfrom
    )
  )
  expect_equal(raster_args$lookat, lookat)
  expect_equal(raster_args$lookfrom, lookfrom)
  expect_equal(raster_args$width, 64L)
  expect_equal(raster_args$height, 65L)
})

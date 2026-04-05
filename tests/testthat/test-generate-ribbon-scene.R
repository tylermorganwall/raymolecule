make_two_chain_model = function() {
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

  expect_true(all(vapply(scene$materials, function(x) x[[1]]$diffuse_texname, character(1)) == texture))
  expect_true(length(scene$texcoords) > 0)
})

test_that("generate_ribbon_scene returns valid rayvertex and rayrender scenes", {
  model = make_two_chain_model()

  raster_scene = generate_ribbon_scene(model, pathtrace = FALSE, subdivisions = 2)
  path_scene = generate_ribbon_scene(model, pathtrace = TRUE, subdivisions = 2)

  expect_s3_class(raster_scene, "ray_mesh")
  expect_s3_class(path_scene, "ray_scene")
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

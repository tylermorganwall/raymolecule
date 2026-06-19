test_that("atom, bond, and full scene generators return raymesh scenes", {
  model = read_sdf(get_example_molecule("caffeine"))

  expect_s3_class(generate_atom_scene(model), "ray_mesh")
  expect_s3_class(generate_bond_scene(model), "ray_mesh")
  expect_s3_class(generate_full_scene(model), "ray_mesh")
})

test_that("pathtrace is selected by render_model rather than scene generators", {
  expect_false("pathtrace" %in% names(formals(generate_atom_scene)))
  expect_false("pathtrace" %in% names(formals(generate_bond_scene)))
  expect_false("pathtrace" %in% names(formals(generate_full_scene)))
  expect_false("pathtrace" %in% names(formals(generate_ribbon_scene)))
  expect_true("pathtrace" %in% names(formals(render_model)))
})

test_that("scene generators accept rayrender material args", {
  expect_true("material_args" %in% names(formals(generate_atom_scene)))
  expect_true("material_args" %in% names(formals(generate_bond_scene)))
  expect_true("material_args" %in% names(formals(generate_full_scene)))
  expect_true("material_args" %in% names(formals(generate_ribbon_scene)))
})

test_that("material_args passes supported material arguments", {
  mesh_material = rayrender_material_to_vertex_material(
    rayrender::dielectric,
    material_args = list(
      refraction = 1.8,
      bump_texture = "height.png",
      bump_intensity = 2
    )
  )

  expect_equal(mesh_material$ior, 1.8)
  expect_equal(mesh_material$bump_texname, "height.png")
  expect_equal(mesh_material$bump_intensity, 2)
})

test_that("glossy material args affect converted mesh material", {
  matte_material = rayrender_material_to_vertex_material(
    rayrender::glossy,
    material_args = list(gloss = 0, reflectance = 0)
  )
  shiny_material = rayrender_material_to_vertex_material(
    rayrender::glossy,
    material_args = list(gloss = 1, reflectance = 0.5)
  )

  expect_equal(matte_material$reflection_intensity, 0)
  expect_equal(matte_material$specular_intensity, 0)
  expect_equal(matte_material$reflection_sharpness, 0)
  expect_equal(shiny_material$reflection_intensity, 0.5)
  expect_equal(shiny_material$specular_intensity, 0.5)
  expect_equal(shiny_material$reflection_sharpness, 1)
})

test_that("pathtraced raymeshes override attached rayrender materials", {
  model = read_sdf(get_example_molecule("benzene"))
  scene = generate_full_scene(
    model,
    material = rayrender::glossy,
    material_args = list(gloss = 0, reflectance = 0)
  )

  rayrender_scene = raymesh_scene_to_rayrender_scene(scene)

  override_materials = vapply(
    rayrender_scene$shape_info,
    function(shape_info) shape_info$shape_properties$override_material,
    logical(1)
  )
  reflectance = vapply(
    rayrender_scene$material,
    function(material) material$glossyinfo[[1]][4],
    numeric(1)
  )

  expect_true(all(override_materials))
  expect_true(all(reflectance == 0))
  expect_lt(nrow(rayrender_scene), length(scene$shapes))
})

test_that("pathtraced raymeshes fall back without attached rayrender materials", {
  model = read_sdf(get_example_molecule("benzene"))
  scene = generate_full_scene(
    model,
    material_vertex = rayvertex::material_list(type = "phong")
  )

  rayrender_scene = raymesh_scene_to_rayrender_scene(scene)

  override_materials = vapply(
    rayrender_scene$shape_info,
    function(shape_info) shape_info$shape_properties$override_material,
    logical(1)
  )
  expect_false(any(override_materials))
})

test_that("package material args override matching material_args", {
  mesh_material = rayrender_material_to_vertex_material(
    rayrender::diffuse,
    material_args = list(
      color = "red",
      image_texture = "texture.png",
      sigma = 0.4
    ),
    package_args = list(
      color = "white",
      image_texture = "package.png"
    )
  )

  expect_equal(mesh_material$diffuse, c(1, 1, 1))
  expect_equal(mesh_material$diffuse_texname, "package.png")
  expect_equal(mesh_material$type, "diffuse")
  expect_equal(mesh_material$sigma, rayrender::diffuse(sigma = 0.4)[[1]]$sigma)
})

test_that("material_args validates material arguments", {
  expect_error(
    rayrender_material_to_vertex_material(
      rayrender::glossy,
      material_args = c(gloss = 0.5)
    ),
    "must be a named list"
  )
  expect_error(
    rayrender_material_to_vertex_material(
      rayrender::glossy,
      material_args = list(0.5)
    ),
    "must be a named list"
  )
  expect_error(
    rayrender_material_to_vertex_material(
      rayrender::glossy,
      material_args = list(glossiness = 0.5)
    ),
    "Unknown rayrender material argument"
  )
})

test_that("current rayrender diffuse materials resolve to rayvertex materials", {
  mesh_material = rayrender_material_to_vertex_material(rayrender::diffuse)

  expect_equal(mesh_material$type, "diffuse")
  expect_equal(mesh_material$specular_intensity, 0)
})

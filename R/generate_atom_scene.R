#' Build Scene (atoms only)
#'
#' Builds a scene containing only atom spheres from a PDB or SDF model.
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the inter-atom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param material Default `rayrender::glossy`. Optional rayrender material
#'   used to initialize the mesh material when `material_vertex` is not
#'   supplied. Must be either `glossy`, `diffuse`, or `dielectric`.
#' @param material_args Default `list()`. Named list of additional arguments
#'   passed to `material`. Arguments supplied by raymolecule for colors and
#'   textures override entries with the same names. For example, use
#'   `list(gloss = 0.35, reflectance = 0.12)` with `rayrender::glossy`, or
#'   `list(sigma = 0.4)` with `rayrender::diffuse`.
#' @param material_vertex Default `rayvertex::material_list()`. Mesh material.
#' `diffuse`/`ambient` colors and `ambient_intensity` are determined automatically, but all other material
#' properties can be changed.
#'
#' @return Raymesh scene containing only the atoms in a molecule/protein.
#' @importFrom rayrender glossy
#' @importFrom rayvertex material_list sphere_mesh scene_from_list
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' atom_model = read_sdf(get_example_molecule("benzene"))
#'
#' # Start with a centered raster atom scene using the default atom colors.
#' atom_model |>
#'   generate_atom_scene() |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey15"
#'   )
#'
#' # This version shifts and scales the same atoms, then uses a diffuse mesh
#' # material for the spheres before pathtracing.
#' atom_model |>
#'   generate_atom_scene(
#'     x = -1,
#'     y = 0,
#'     z = 1,
#'     scale = 0.75,
#'     center = TRUE,
#'     material = rayrender::diffuse,
#'     material_args = list(sigma = 0.4)
#'   ) |>
#'   render_model(pathtrace = TRUE, width = 800, height = 800, samples = 32)
#'
#' # The toon material changes only the raster shader. Keeping `center = TRUE`
#' # makes the result easy to compare to the baseline scene above.
#' toon_material = rayvertex::material_list(
#'   type = "toon_phong",
#'   toon_levels = 3,
#'   toon_outline_width = 0.08
#' )
#' caffeine_model = read_sdf(get_example_molecule("caffeine"))
#' caffeine_model |>
#'   generate_atom_scene(
#'     center = TRUE,
#'     material_vertex = toon_material
#'   ) |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey15"
#'   )
generate_atom_scene = function(
  model,
  x = 0,
  y = 0,
  z = 0,
  scale = 1,
  center = TRUE,
  material = rayrender::glossy,
  material_args = list(),
  material_vertex = material_list(type = "phong")
) {
  material_template = resolve_mesh_material(
    material = material,
    material_vertex = material_vertex,
    use_material = !missing(material) && missing(material_vertex),
    material_args = material_args,
    package_args = list(color = "white", image_texture = "")
  )
  atoms = model$atoms
  atoms$x = atoms$x * scale
  atoms$y = atoms$y * scale
  atoms$z = atoms$z * scale
  if (center) {
    atoms$x = atoms$x - mean(range(atoms$x))
    atoms$y = atoms$y - mean(range(atoms$y))
    atoms$z = atoms$z - mean(range(atoms$z))
  }
  atoms$x = atoms$x + x
  atoms$y = atoms$y + y
  atoms$z = atoms$z + z

  scene = list()
  for (i in seq_len(nrow(atoms))) {
    if (atoms$type[i] != "C") {
      atomcol = PeriodicTable::atomColor(atoms$type[i])
    } else {
      atomcol = "grey5"
    }

    material_atom = material_template
    material_atom$diffuse = as.numeric(grDevices::col2rgb(atomcol)) / 255
    material_atom$ambient = as.numeric(grDevices::col2rgb(atomcol)) / 255
    material_atom$ambient_intensity = 0.3
    material_atom = update_rayrender_material_package_args(
      material_atom,
      package_args = list(color = atomcol, image_texture = "")
    )
    atomsize = (PeriodicTable::mass(atoms$type[i]) / 14)^(1 / 3)
    scene[[i]] = sphere_mesh(
      position = c(atoms$x[i], atoms$y[i], atoms$z[i]),
      radius = atomsize / 2,
      low_poly = FALSE,
      material = material_atom
    )
  }
  return(rayvertex::scene_from_list(scene))
}

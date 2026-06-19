#' Build  Scene (bonds + atoms)
#'
#' Builds a combined atom and bond scene from a PDB or SDF model.
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the interatom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param force_single_bonds Default `FALSE`. Whether to force all bonds to show as a single connection.
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
#'
#' @return Raymesh scene
#' @importFrom rayrender glossy
#' @importFrom rayvertex material_list sphere_mesh segment_mesh scene_from_list
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' molecule_model = read_sdf(get_example_molecule("caffeine"))
#'
#' # Start with a centered raster scene that combines atom spheres and bond
#' # geometry using the default atom colors.
#' molecule_model |>
#'   generate_full_scene(force_single_bonds = TRUE) |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey12"
#'   )
#'
#' # This version changes the scale, keeps the model centered, and uses a
#' # diffuse mesh material for both atoms and bonds before pathtracing.
#' molecule_model |>
#'   generate_full_scene(
#'     x = 0,
#'     y = 0,
#'     z = 0,
#'     scale = 0.75,
#'     center = TRUE,
#'     force_single_bonds = TRUE,
#'     material = rayrender::diffuse,
#'     material_args = list(sigma = 0.3)
#'   ) |>
#'   render_model(pathtrace = TRUE, width = 800, height = 800, samples = 32)
#'
#' # A toon material changes the raster shader and pass FSAA to render_model()
#' # so rayvertex adds anti-aliasing.
#' shiny_toon_material = rayvertex::material_list(
#'   type = "toon_phong",
#'   toon_levels = 3,
#'   toon_outline_width = 10,
#'   toon_outline_color = "white"
#' )
#' morphine_model = read_sdf(get_example_molecule("morphine"))
#' morphine_model |>
#'   generate_full_scene(
#'     material_vertex = shiny_toon_material
#'   ) |>
#'   render_model(
#'     fsaa = 2,
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey50"
#'   )
generate_full_scene = function(
  model,
  x = 0,
  y = 0,
  z = 0,
  scale = 1,
  center = TRUE,
  force_single_bonds = FALSE,
  material = rayrender::glossy,
  material_args = list(),
  material_vertex = material_list(type = "phong")
) {
  material_vertex = resolve_mesh_material(
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
  cntr = 1
  for (i in seq_len(nrow(atoms))) {
    if (atoms$type[i] != "C") {
      atomcol = PeriodicTable::atomColor(atoms$type[i])
    } else {
      atomcol = "grey5"
    }
    material_atom = material_vertex
    material_atom$diffuse = convert_color(atomcol)
    material_atom$ambient = convert_color(atomcol)
    material_atom$ambient_intensity = 0.3
    material_atom = update_rayrender_material_package_args(
      material_atom,
      package_args = list(color = atomcol, image_texture = "")
    )
    atomsize = (PeriodicTable::mass(atoms$type[i]) / 14)^(1 / 3)
    scene[[cntr]] = sphere_mesh(
      position = c(atoms$x[i], atoms$y[i], atoms$z[i]),
      radius = atomsize / 2,
      low_poly = FALSE,
      material = material_atom
    )
    cntr = cntr + 1
  }

  get_segment_offsets = function(start, end, bond_order) {
    zero_offset = c(0, 0, 0)
    if (bond_order == 1 || force_single_bonds) {
      return(list(zero_offset))
    }
    dir = end - start
    onb = onb_from_w(dir)
    if (bond_order == 2) {
      offset = onb[3, ] / 8
      return(list(offset, -offset))
    }
    if (bond_order == 3) {
      offset = onb[3, ] / 4
      return(list(offset, -offset, zero_offset))
    }
    return(NULL)
  }

  material_bond = material_vertex
  material_bond$diffuse = convert_color("grey33")
  material_bond$ambient = convert_color("grey33")
  material_bond$ambient_intensity = 0.3
  material_bond = update_rayrender_material_package_args(
    material_bond,
    package_args = list(color = "grey33", image_texture = "")
  )
  bonds = model$bonds
  for (i in seq_len(nrow(bonds))) {
    bond1 = atoms$index == bonds[i, 1]
    bond2 = atoms$index == bonds[i, 2]
    if (!any(bond1) || !any(bond2)) {
      next
    }

    start = as.numeric(atoms[bond1, c("x", "y", "z")])
    end = as.numeric(atoms[bond2, c("x", "y", "z")])
    segment_offsets = get_segment_offsets(start, end, bonds[i, 3])
    if (is.null(segment_offsets)) {
      next
    }

    for (offset in segment_offsets) {
      scene[[cntr]] = segment_mesh(
        start = start + offset,
        end = end + offset,
        radius = 1 / 10,
        material = material_bond
      )
      cntr = cntr + 1
    }
  }

  return(rayvertex::scene_from_list(scene))
}

#' Build Scene (bonds only)
#'
#' Builds a scene containing only bond geometry from a PDB or SDF model.
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
#' @param material_vertex Default `material_list(diffuse="grey33",ambient="grey33",type="phong", ambient_intensity=0.3)`.
#' Mesh material to use for the bonds.
#'
#' @return Raymesh scene containing only the connections between atoms in a molecule/protein.
#' @importFrom rayrender glossy
#' @importFrom rayvertex material_list sphere_mesh segment_mesh scene_from_list
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' bond_model = read_sdf(get_example_molecule("benzene"))
#'
#' # Start with a centered bond scene using the molecule's recorded bond
#' # orders and default bond material.
#' bond_model |>
#'   generate_bond_scene() |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey10"
#'   )
#'
#' # This version reduces the spacing, forces every connection to a single
#' # bond, and lights the ligand from above and below while pathtracing.
#' bond_model |>
#'   generate_bond_scene(
#'     x = 0,
#'     y = 0,
#'     z = 0,
#'     scale = 0.7,
#'     center = TRUE,
#'     force_single_bonds = TRUE,
#'     material = rayrender::glossy,
#'     material_args = list(gloss = 0.35)
#'   ) |>
#'   render_model(
#'     pathtrace = TRUE,
#'     lights = "both",
#'     width = 800,
#'     height = 800,
#'     samples = 32
#'   )
#'
#' # A custom rayvertex material changes the raster bond color and ambient
#' # contribution while keeping the same geometry.
#' bond_material = rayvertex::material_list(
#'   diffuse = "grey85",
#'   ambient = "grey25",
#'   type = "phong",
#'   ambient_intensity = 0.4
#' )
#' cinnemaldehyde_model = read_sdf(get_example_molecule("cinnemaldehyde"))
#' cinnemaldehyde_model |>
#'   generate_bond_scene(
#'     material_vertex = bond_material
#'   ) |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey10"
#'   )
generate_bond_scene = function(
  model,
  x = 0,
  y = 0,
  z = 0,
  scale = 1,
  center = TRUE,
  force_single_bonds = FALSE,
  material = rayrender::glossy,
  material_args = list(),
  material_vertex = material_list(
    diffuse = "grey33",
    ambient = "grey33",
    type = "phong",
    ambient_intensity = 0.3
  )
) {
  material_vertex = resolve_mesh_material(
    material = material,
    material_vertex = material_vertex,
    use_material = !missing(material) && missing(material_vertex),
    material_args = material_args,
    package_args = list(color = "grey33", image_texture = "")
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

  get_bond_offsets = function(start, end, bond_order) {
    zero_offset = c(0, 0, 0)
    if (bond_order == 1 || force_single_bonds) {
      return(list(
        segment_offsets = list(zero_offset, zero_offset),
        endpoint_offsets = list(zero_offset)
      ))
    }
    dir = end - start
    onb = onb_from_w(dir)
    if (bond_order == 2) {
      offset = onb[3, ] / 8
      return(list(
        segment_offsets = list(offset, -offset),
        endpoint_offsets = list(offset, -offset)
      ))
    }
    if (bond_order == 3) {
      offset = onb[3, ] / 4
      return(list(
        segment_offsets = list(offset, -offset, zero_offset),
        endpoint_offsets = list(offset, -offset, zero_offset)
      ))
    }
    return(NULL)
  }

  bonds = model$bonds
  bond_scene = list()
  counter = 1
  for (i in seq_len(nrow(bonds))) {
    bond1 = atoms$index == bonds[i, 1]
    bond2 = atoms$index == bonds[i, 2]
    if (!any(bond1) || !any(bond2)) {
      next
    }

    start = as.numeric(atoms[bond1, c("x", "y", "z")])
    end = as.numeric(atoms[bond2, c("x", "y", "z")])
    offsets = get_bond_offsets(start, end, bonds[i, 3])
    if (is.null(offsets)) {
      next
    }

    for (offset in offsets$segment_offsets) {
      bond_scene[[counter]] = segment_mesh(
        start = start + offset,
        end = end + offset,
        radius = 1 / 10,
        material = material_vertex
      )
      counter = counter + 1
    }

    for (offset in offsets$endpoint_offsets) {
      bond_scene[[counter]] = sphere_mesh(
        position = start + offset,
        radius = 1 / 10,
        material = material_vertex
      )
      counter = counter + 1
      bond_scene[[counter]] = sphere_mesh(
        position = end + offset,
        radius = 1 / 10,
        material = material_vertex
      )
      counter = counter + 1
    }
  }

  return(rayvertex::scene_from_list(bond_scene))
}

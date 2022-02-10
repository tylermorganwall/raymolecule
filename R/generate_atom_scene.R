#' Build Scene (atoms only)
#'
#' Reads an SDF file and extracts the 3D molecule model
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the inter-atom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param pathtrace Default `TRUE`. If `FALSE`, the `rayvertex` package will be used to render the scene.
#' @param material Default `rayrender::glossy`. Rayrender material to use when `pathtrace = TRUE`. Must be either `glossy`, `diffuse`, or `dielectric`.
#' @param material_vertex Default `rayvertex::material_list()`. Material to use when `pathtrace = FALSE`.
#' `diffuse`/`ambient` colors and `ambient_intensity` are determined automatically, but all other material
#' properties can be changed.
#'
#' @return Rayrender/rayvertex scene containing only the atoms in a molecule/protein.
#' @importFrom rayrender render_scene glossy sphere segment add_object
#' @importFrom rayvertex rasterize_scene material_list sphere_mesh segment_mesh add_shape
#' @export
#'
#' @examples
#' #Generate a scene with caffeine molecule with just the atoms
#'\donttest{
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_atom_scene() %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'
#' #Generate a rayvertex scene, using toon shading
#' shiny_toon_material = rayvertex::material_list(type="toon_phong",
#'                                                toon_levels=3,
#'                                                toon_outline_width=0.1)
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_atom_scene(pathtrace=FALSE, material_vertex = shiny_toon_material) %>%
#'   render_model(background="white")
#'
#'#Generate a scene with caffeine, reducing the inter-atom spacing
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_atom_scene(scale=0.5) %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'}
generate_atom_scene = function(model, x=0, y=0, z=0, scale = 1, center = TRUE,
                               pathtrace = TRUE,
                               material = rayrender::glossy,
                               material_vertex = material_list(type="phong")) {
  mat_info = material()
  if(!mat_info$type %in% c("glossy","diffuse", "dielectric")) {
    stop("material() must be either `glossy`, `diffuse`, or `dielectric`")
  }
  atoms = model$atoms
  atoms$x = atoms$x * scale
  atoms$y = atoms$y * scale
  atoms$z = atoms$z * scale
  if(center) {
    atoms$x = atoms$x - mean(range(atoms$x))
    atoms$y = atoms$y - mean(range(atoms$y))
    atoms$z = atoms$z - mean(range(atoms$z))
  }
  atoms$x = atoms$x + x
  atoms$y = atoms$y + y
  atoms$z = atoms$z + z
  if(pathtrace) {
    scenelist = list()
    counter = 1
    for(i in 1:nrow(atoms)) {
      if(atoms$type[i] != "C") {
        atomcol = PeriodicTable::atomColor(atoms$type[i])
      } else {
        atomcol = "grey5"
      }
      atomsize = (PeriodicTable::mass(atoms$type[i])/14)^(1/3)

      scenelist[[counter]] = sphere(x=atoms$x[i],y=atoms$y[i],z=atoms$z[i],
                                    radius=atomsize/2,
                                    material=material(color=atomcol))
      counter = counter + 1
    }
    return(do.call(rbind, scenelist))
  } else {
    scene = list()
    for (i in 1:nrow(atoms)) {
      if (atoms$type[i] != "C") {
        atomcol = PeriodicTable::atomColor(atoms$type[i])
      }
      else {
        atomcol = "grey5"
      }
      material_atom = material_vertex
      material_atom$diffuse = atomcol
      material_atom$ambient = atomcol
      material_atom$ambient_intensity = 0.3
      atomsize = (PeriodicTable::mass(atoms$type[i])/14)^(1/3)
      scene = add_shape(scene, sphere_mesh(position = c(atoms$x[i],atoms$y[i], atoms$z[i]),
                                           radius = atomsize/2, low_poly = F,
                                           material = material_atom))
    }
    return(scene)
  }
}

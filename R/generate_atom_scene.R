#' Build Rayrender Scene (atoms only)
#'
#' Reads an SDF file and extracts the 3D molecule model
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the inter-atom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#'
#' @return Rayrender scene (tibble) containing only the atoms in a molecule/protein.
#' @import rayrender
#' @export
#'
#' @examples
#' #Generate a scene with caffeine molecule with just the atoms
#'\donttest{
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_atom_scene() %>%
#'   render_model()
#'
#'#Generate a scene with caffeine, reducing the inter-atom spacing
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_atom_scene(scale=0.5) %>%
#'   render_model()
#'}
generate_atom_scene = function(model, x=0, y=0, z=0, scale = 1, center = TRUE) {
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
                                  material=glossy(color=atomcol))
    counter = counter + 1
  }
  return(do.call(rbind, scenelist))
}

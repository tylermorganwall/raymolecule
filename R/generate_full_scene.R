#' Build Rayrender Scene (bonds + atoms)
#'
#' Reads an SDF file and extracts the 3D molecule model
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param scale Default `1`. Amount to scale the interatom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param force_single_bonds Default `FALSE`. Whether to force all bonds to show as a single connection.
#'
#' @return Rayrender scene (tibble) containing the atoms and bonds in a molecule/protein.
#' @import rayrender
#' @export
#'
#' @examples
#' # Generate a scene with caffeine molecule
#'\donttest{
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' # Generate a scene with morphine, increasing the inter-atom spacing
#' get_example_molecule("tubocurarine_chloride") %>%
#'   read_sdf() %>%
#'   generate_full_scene(scale=1.5) %>%
#'   render_model()
#'
#' # Force bonds to appear as a single link (to focus purely on the shape of the molecule)
#' get_example_molecule("tubocurarine_chloride") %>%
#'   read_sdf() %>%
#'   generate_full_scene(force_single_bonds = TRUE) %>%
#'   render_model()
#'}
generate_full_scene = function(model, scale = 1, center = TRUE, force_single_bonds = FALSE) {
  atoms = model$atoms
  atoms$x = atoms$x * scale
  atoms$y = atoms$y * scale
  atoms$z = atoms$z * scale
  if(center) {
    atoms$x = atoms$x - mean(range(atoms$x))
    atoms$y = atoms$y - mean(range(atoms$y))
    atoms$z = atoms$z - mean(range(atoms$z))
  }
  bonds = model$bonds
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
  bondlist = list()
  counter = 1
  for(i in 1:nrow(bonds)) {
    bond1 = atoms$index == bonds[i,1]
    bond2 = atoms$index == bonds[i,2]
    if(bonds[i,3] == 1 || force_single_bonds) {
      if(any(bond1) && any(bond2)) {
        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3]),
                                      end = as.numeric(atoms[bond2,1:3]),
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1
      }

    } else if(bonds[i,3] == 2) {
      if(any(bond1) && any(bond2)) {
        dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
        onb = onb_from_w(dir)
        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])+onb[3,]/8,
                                      end = as.numeric(atoms[bond2,1:3])+onb[3,]/8,
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1

        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])-onb[3,]/8,
                                      end = as.numeric(atoms[bond2,1:3])-onb[3,]/8,
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1
      }
    } else if(bonds[i,3] == 3) {
      if(any(bond1) && any(bond2)) {
        dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
        onb = onb_from_w(dir)
        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])+onb[3,]/4,
                                      end = as.numeric(atoms[bond2,1:3])+onb[3,]/4,
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1

        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])-onb[3,]/4,
                                      end = as.numeric(atoms[bond2,1:3])-onb[3,]/4,
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1
        bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3]),
                                      end = as.numeric(atoms[bond2,1:3]),
                                      radius=1/10,
                                      material=glossy(color="grey33"))
        counter = counter + 1
      }
    }
  }

  atom_scene = do.call(rbind, scenelist)
  bond_scene = do.call(rbind, bondlist)
  return(add_object(atom_scene,bond_scene))
}

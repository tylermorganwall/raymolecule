#' Build Scene (bonds only)
#'
#' Reads an SDF file and extracts the 3D molecule model
#'
#' @param model Model extracted from a PDB or SDF file.
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the interatom spacing.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param force_single_bonds Default `FALSE`. Whether to force all bonds to show as a single connection.
#' @param pathtrace Default `TRUE`. If `FALSE`, the `rayvertex` package will be used to render the scene.
#' @param material Default `rayrender::glossy`. Rayrender material to use when `pathtrace = TRUE`. Must be either `glossy`, `diffuse`, or `dielectric`.
#' @param material_vertex Default `material_list(diffuse="grey33",ambient="grey33",type="phong", ambient_intensity=0.3)`.
#' Material to use for the bonds when `pathtrace = FALSE`.
#'
#' @return Rayrender/rayvertex scene containing only the connections between atoms in a molecule/protein.
#' @importFrom rayrender render_scene glossy sphere segment add_object
#' @importFrom rayvertex rasterize_scene material_list sphere_mesh segment_mesh add_shape
#' @export
#'
#' @examples
#' #Generate a scene with benzene molecule with just the atoms
#'\donttest{
#' get_example_molecule("benzene") %>%
#'   read_sdf() %>%
#'   generate_bond_scene() %>%
#'   render_model(lights = "both", samples=256,sample_method="sobol_blue")
#'
#'#Force single bonds to just show the shape of the molecule
#' get_example_molecule("benzene") %>%
#'   read_sdf() %>%
#'   generate_bond_scene(force_single_bonds = TRUE) %>%
#'   render_model(lights = "both", samples=256,sample_method="sobol_blue")
#'
#'#Generate a scene with PFOA, reducing the inter-atom spacing
#' get_example_molecule("pfoa") %>%
#'   read_sdf() %>%
#'   generate_bond_scene(scale=0.3,force_single_bonds = TRUE) %>%
#'   render_model(lights = "both", samples=256,sample_method="sobol_blue")
#'}
generate_bond_scene = function(model, x=0, y=0, z=0, scale = 1, center = TRUE,
                               force_single_bonds = FALSE, pathtrace = TRUE,
                               material = rayrender::glossy,
                               material_vertex = material_list(diffuse="grey33",
                                                               ambient="grey33",
                                                               type = "phong",
                                                               ambient_intensity=0.3)) {
  mat_info = material()
  mat_info = mat_info[[1]]
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

  bonds = model$bonds
  if(pathtrace) {
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
                                        material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3]),
                                        end = as.numeric(atoms[bond2,1:3]),
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1

          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1]),
                                       y = as.numeric(atoms[bond1,2]),
                                       z = as.numeric(atoms[bond1,3]),
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1]),
                                       y = as.numeric(atoms[bond2,2]),
                                       z = as.numeric(atoms[bond2,3]),
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
        }
      } else if(bonds[i,3] == 2) {
        if(any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
          onb = onb_from_w(dir)
          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])+onb[3,]/8,
                                        end = as.numeric(atoms[bond2,1:3])+onb[3,]/8,
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1


          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])-onb[3,]/8,
                                        end = as.numeric(atoms[bond2,1:3])-onb[3,]/8,
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1])+onb[3,1]/8,
                                       y = as.numeric(atoms[bond1,2])+onb[3,2]/8,
                                       z = as.numeric(atoms[bond1,3])+onb[3,3]/8,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1])+onb[3,1]/8,
                                       y = as.numeric(atoms[bond2,2])+onb[3,2]/8,
                                       z = as.numeric(atoms[bond2,3])+onb[3,3]/8,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1])-onb[3,1]/8,
                                       y = as.numeric(atoms[bond1,2])-onb[3,2]/8,
                                       z = as.numeric(atoms[bond1,3])-onb[3,3]/8,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1])-onb[3,1]/8,
                                       y = as.numeric(atoms[bond2,2])-onb[3,2]/8,
                                       z = as.numeric(atoms[bond2,3])-onb[3,3]/8,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
        }
      } else if(bonds[i,3] == 3) {
        if(any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
          onb = onb_from_w(dir)
          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])+onb[3,]/4,
                                        end = as.numeric(atoms[bond2,1:3])+onb[3,]/4,
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1

          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3])-onb[3,]/4,
                                        end = as.numeric(atoms[bond2,1:3])-onb[3,]/4,
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = segment(start = as.numeric(atoms[bond1,1:3]),
                                        end = as.numeric(atoms[bond2,1:3]),
                                        radius=1/10,
                                        material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1])+onb[3,1]/4,
                                       y = as.numeric(atoms[bond1,2])+onb[3,2]/4,
                                       z = as.numeric(atoms[bond1,3])+onb[3,3]/4,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1])+onb[3,1]/4,
                                       y = as.numeric(atoms[bond2,2])+onb[3,2]/4,
                                       z = as.numeric(atoms[bond2,3])+onb[3,3]/4,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1])-onb[3,1]/4,
                                       y = as.numeric(atoms[bond1,2])-onb[3,2]/4,
                                       z = as.numeric(atoms[bond1,3])-onb[3,3]/4,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1])-onb[3,1]/4,
                                       y = as.numeric(atoms[bond2,2])-onb[3,2]/4,
                                       z = as.numeric(atoms[bond2,3])-onb[3,3]/4,
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond1,1]),
                                       y = as.numeric(atoms[bond1,2]),
                                       z = as.numeric(atoms[bond1,3]),
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
          bondlist[[counter]] = sphere(x = as.numeric(atoms[bond2,1]),
                                       y = as.numeric(atoms[bond2,2]),
                                       z = as.numeric(atoms[bond2,3]),
                                       radius=1/10,
                                       material=material(color="grey33"))
          counter = counter + 1
        }
      }
    }

    bond_scene = do.call(rbind, bondlist)
  } else {
    bond_scene = list()
    counter = 1

    for(i in 1:nrow(bonds)) {
      bond1 = atoms$index == bonds[i,1]
      bond2 = atoms$index == bonds[i,2]
      if(bonds[i,3] == 1 || force_single_bonds) {
        if(any(bond1) && any(bond2)) {
          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3]),
                                        end = as.numeric(atoms[bond2,1:3]),
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3]),
                                        end = as.numeric(atoms[bond2,1:3]),
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1

          bond_scene = add_shape(bond_scene, sphere_mesh(position=c(as.numeric(atoms[bond1,1]),
                                       as.numeric(atoms[bond1,2]),
                                       as.numeric(atoms[bond1,3])),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1]),
                                       as.numeric(atoms[bond2,2]),
                                       as.numeric(atoms[bond2,3])),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
        }
      } else if(bonds[i,3] == 2) {
        if(any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
          onb = onb_from_w(dir)
          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3])+onb[3,]/8,
                                        end = as.numeric(atoms[bond2,1:3])+onb[3,]/8,
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1


          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3])-onb[3,]/8,
                                        end = as.numeric(atoms[bond2,1:3])-onb[3,]/8,
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond1,1])+onb[3,1]/8,
                                       as.numeric(atoms[bond1,2])+onb[3,2]/8,
                                       as.numeric(atoms[bond1,3])+onb[3,3]/8),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1])+onb[3,1]/8,
                                       as.numeric(atoms[bond2,2])+onb[3,2]/8,
                                       as.numeric(atoms[bond2,3])+onb[3,3]/8),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond1,1])-onb[3,1]/8,
                                       as.numeric(atoms[bond1,2])-onb[3,2]/8,
                                       as.numeric(atoms[bond1,3])-onb[3,3]/8),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1])-onb[3,1]/8,
                                       as.numeric(atoms[bond2,2])-onb[3,2]/8,
                                       as.numeric(atoms[bond2,3])-onb[3,3]/8),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
        }
      } else if(bonds[i,3] == 3) {
        if(any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2,1:3])-as.numeric(atoms[bond1,1:3])
          onb = onb_from_w(dir)
          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3])+onb[3,]/4,
                                        end = as.numeric(atoms[bond2,1:3])+onb[3,]/4,
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1

          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3])-onb[3,]/4,
                                        end = as.numeric(atoms[bond2,1:3])-onb[3,]/4,
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, segment_mesh(start = as.numeric(atoms[bond1,1:3]),
                                        end = as.numeric(atoms[bond2,1:3]),
                                        radius=1/10,
                                        material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond1,1])+onb[3,1]/4,
                                       as.numeric(atoms[bond1,2])+onb[3,2]/4,
                                       as.numeric(atoms[bond1,3])+onb[3,3]/4),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1])+onb[3,1]/4,
                                       as.numeric(atoms[bond2,2])+onb[3,2]/4,
                                       as.numeric(atoms[bond2,3])+onb[3,3]/4),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond1,1])-onb[3,1]/4,
                                       as.numeric(atoms[bond1,2])-onb[3,2]/4,
                                       as.numeric(atoms[bond1,3])-onb[3,3]/4),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1])-onb[3,1]/4,
                                       as.numeric(atoms[bond2,2])-onb[3,2]/4,
                                       as.numeric(atoms[bond2,3])-onb[3,3]/4),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond1,1]),
                                       as.numeric(atoms[bond1,2]),
                                       as.numeric(atoms[bond1,3])),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
          bond_scene = add_shape(bond_scene, sphere_mesh(position = c(as.numeric(atoms[bond2,1]),
                                       as.numeric(atoms[bond2,2]),
                                       as.numeric(atoms[bond2,3])),
                                       radius=1/10,
                                       material=material_vertex))
          counter = counter + 1
        }
      }
    }
  }
  return(bond_scene)
}

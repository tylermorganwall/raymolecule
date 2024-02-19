#' Build  Scene (bonds + atoms)
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
#' @param material_vertex Default `rayvertex::material_list()`. Material to use when `pathtrace = FALSE`.
#' `diffuse`/`ambient` colors and `ambient_intensity` are determined automatically, but all other material
#' properties can be changed.
#'
#'
#' @return Rayrender/rayvertex scene
#' @importFrom rayrender render_scene glossy sphere segment add_object
#' @importFrom rayvertex rasterize_scene material_list sphere_mesh segment_mesh add_shape
#' @export
#'
#' @examples
#' # Generate a scene with caffeine molecule
#'\donttest{
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'
#' #Generate a rayvertex scene with a custom material
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene(pathtrace=FALSE, material_vertex=rayvertex::material_list(type="phong")) %>%
#'   render_model(background="grey33")
#'
#' #Generate a rayvertex scene, using toon shading
#' shiny_toon_material = rayvertex::material_list(type="toon_phong",
#'                                                toon_levels=3,
#'                                                toon_outline_width=0.1)
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene(pathtrace=FALSE, material_vertex=shiny_toon_material) %>%
#'   render_model(background="grey66")
#'
#' # Generate a scene with morphine, increasing the inter-atom spacing
#' get_example_molecule("tubocurarine_chloride") %>%
#'   read_sdf() %>%
#'   generate_full_scene(scale=1.5) %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'
#'
#' # Force bonds to appear as a single link (to focus purely on the shape of the molecule)
#' get_example_molecule("tubocurarine_chloride") %>%
#'   read_sdf() %>%
#'   generate_full_scene(force_single_bonds = TRUE) %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'
#'}
generate_full_scene = function(model, x=0,y=0,z=0, scale = 1, center = TRUE, pathtrace = TRUE,
                               force_single_bonds = FALSE, material = rayrender::glossy,
                               material_vertex = material_list(type="phong")) {
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
        }
      }
    }

    atom_scene = do.call(rbind, scenelist)
    bond_scene = do.call(rbind, bondlist)
    return(add_object(atom_scene,bond_scene))
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
      material_atom$diffuse = convert_color(atomcol)
      material_atom$ambient = convert_color(atomcol)
      material_atom$ambient_intensity = 0.3
      atomsize = (PeriodicTable::mass(atoms$type[i])/14)^(1/3)
      scene = add_shape(scene, sphere_mesh(position = c(atoms$x[i],atoms$y[i], atoms$z[i]),
                                           radius = atomsize/2, low_poly = F,
                                           material = material_atom))
    }
    material_bond = material_vertex
    material_atom$diffuse = convert_color("grey33")
    material_atom$ambient = convert_color("grey33")
    material_atom$ambient_intensity = 0.3
    for (i in 1:nrow(bonds)) {
      bond1 = atoms$index == bonds[i, 1]
      bond2 = atoms$index == bonds[i, 2]
      if (bonds[i, 3] == 1 || force_single_bonds) {
        if (any(bond1) && any(bond2)) {
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1,1:3]), end = as.numeric(atoms[bond2, 1:3]),
                                                radius = 1/10,
                                                material=material_bond))
        }
      }
      else if (bonds[i, 3] == 2) {
        if (any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2, 1:3]) - as.numeric(atoms[bond1, 1:3])
          onb = onb_from_w(dir)
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1, 1:3]) + onb[3, ]/8, end = as.numeric(atoms[bond2, 1:3]) + onb[3, ]/8,
                                                radius = 1/10,
                                                material = material_bond))
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1, 1:3]) - onb[3, ]/8, end = as.numeric(atoms[bond2, 1:3]) - onb[3, ]/8,
                                                radius = 1/10,
                                                material = material_bond))
        }
      }
      else if (bonds[i, 3] == 3) {
        if (any(bond1) && any(bond2)) {
          dir = as.numeric(atoms[bond2, 1:3]) - as.numeric(atoms[bond1,1:3])
          onb = onb_from_w(dir)
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1,1:3]) + onb[3, ]/4, end = as.numeric(atoms[bond2, 1:3]) + onb[3, ]/4,
                                                radius = 1/10,
                                                material = material_bond))
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1, 1:3]) - onb[3, ]/4, end = as.numeric(atoms[bond2, 1:3]) - onb[3, ]/4,
                                                radius = 1/10,
                                                material = material_bond))
          scene = add_shape(scene, segment_mesh(start = as.numeric(atoms[bond1, 1:3]), end = as.numeric(atoms[bond2, 1:3]),
                                                radius = 1/10,
                                                material = material_bond))
        }
      }
    }
  }
  return(scene)
}


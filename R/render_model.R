#' Render Molecule Model
#'
#' Automatically plots the molecule with a camera position and field of view that includes the full model.
#' For more control over the scene, pass the scene to `rayrender::render_scene()` or `rayvertex::rasterize_scene()` and specify
#' the camera position manually. Note: spheres and cylinders in the scene are used to automatically
#' compute the field of view of the scene--if rendering with rayrender, adding additional sphere (e.g. with `rayrender::generate_ground()`)
#' will change this calculation. Use `rayrender::render_scene()` instead if this is a problem.
#'
#' @param scene `rayrender` scene of molecule model.
#' @param fov Default `NULL`, automatically calculated. Camera field of view.
#' @param angle Default `c(0,0,0)`. Degrees to rotate the model around the X, Y, and Z axes. If this
#' is a single number, it will be taken as the Y axis rotation.
#' @param order_rotation Default `c(1,2,3)`. What order to apply the rotations specified in `angle`.
#' @param lights Default `top`. If `none`, removes all lights. If `bottom`, lights scene with light
#' underneath model. If `both`, adds lights both above and below model. This can also be a matrix of
#' light information generated with `rayvertex`.
#' @param lightintensity Default `80`. Light intensity for pathtraced scenes.
#' @param ... Other arguments to pass to `rayrender::render_scene()` or `rayvertex::rasterize_scene()`
#'
#' @return Rendered image
#' @importFrom rayrender render_scene glossy sphere segment add_object light group_objects
#' @importFrom rayvertex rasterize_scene material_list sphere_mesh segment_mesh add_shape directional_light add_light rotate_mesh
#' @export
#'
#' @examples
#' # Generate a scene with caffeine molecule with just the atoms
#'\donttest{
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model(samples=256,sample_method="sobol_blue")
#'
#' #Light the example from below as well
#' get_example_molecule("caffeine") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model(lights = "both", samples=256,sample_method="sobol_blue")
#'
#' #Generate a scene with penicillin, increasing the number of samples and the width/height
#' #for a higher quality render.
#' get_example_molecule("penicillin") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model(lights = "both", samples=256, width=800, height=800,sample_method="sobol_blue")
#'
#' #Render the scene with rayvertex and custom lights
#' get_example_molecule("penicillin") %>%
#'   read_sdf() %>%
#'   generate_full_scene(pathtrace=FALSE) %>%
#'   render_model(width=800, height=800,background="grey66",
#'                lights = rayvertex::directional_light(c(0.2,1,1)))
#'
#' #Rotate the molecule 30 degrees around the y axis, and the 30 degrees around the z axis
#' get_example_molecule("penicillin") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   render_model(lights = "both", samples=256, width=800, height=800,
#'                angle=c(0,30,30),sample_method="sobol_blue")
#'
#'
#'
#' #Add a checkered plane underneath, using rayrender::add_object and rayrender::xz_rect().
#' #We also pass a value to `clamp_value` to minimize fireflies (bright spots).
#' library(rayrender)
#' get_example_molecule("skatole") %>%
#'   read_sdf() %>%
#'   generate_full_scene() %>%
#'   add_object(xz_rect(xwidth=1000,zwidth=1000,y=-4,
#'                      material=diffuse(color="#330000",checkercolor="#770000"))) %>%
#'   render_model(samples=256, width=800, height=800, clamp_value=10,
#'                sample_method="sobol_blue")
#'}
render_model = function(scene, fov = NULL, angle = c(0,0,0), order_rotation = c(1,2,3),
                        lights = "top", lightintensity = 80, ...) {
  if(length(angle) == 1) {
    angle = c(0,angle,0)
  }
  pathtraced = suppressWarnings(is.null(scene$vertices))
  if(pathtraced) {
    scene_model = scene[is.na(scene$lightintensity) &
                        (scene$shape == "cylinder" | scene$shape == "sphere"),]
    bbox_x = range(scene_model$x,na.rm=TRUE)
    bbox_y = range(scene_model$y,na.rm=TRUE)
    bbox_z = range(scene_model$z,na.rm=TRUE)
    spheresizes = scene[(scene$shape == "sphere" & scene$type != "light"),4]
    if(length(spheresizes) > 0) {
      max_sphere_radii = max(spheresizes,na.rm=TRUE)
    } else {
      max_sphere_radii = 0.5
    }

    widest = max(c(abs(bbox_x),abs(bbox_y),abs(bbox_z)))
    offset_dist = widest + widest/5 + max_sphere_radii
  } else {
    bbox_x = range(scene$vertices[,1],na.rm=TRUE)
    bbox_y = range(scene$vertices[,2],na.rm=TRUE)
    bbox_z = range(scene$vertices[,3],na.rm=TRUE)
    widest = max(c(abs(bbox_x),abs(bbox_y),abs(bbox_z)))
    max_sphere_radii = 0
  }
  if(is.null(fov)) {
    fov = atan2(widest+widest/5 + max_sphere_radii, widest*5)/pi*180*2
  }
  if(pathtraced) {
    if(any(angle != 0)) {
      scene = group_objects(scene, angle = angle, order_rotation = order_rotation)
    }
    if(lights != "none") {
      if (lights == "top") {
        light = sphere(x=offset_dist*2,y=offset_dist*2,z=offset_dist*2,
                          radius = widest/2,
                          material = light(intensity=lightintensity)) %>%
          add_object(sphere(x=-offset_dist*2,y=offset_dist*2,z=-offset_dist*2,
                            radius = widest/2,
                            material = light(intensity=lightintensity)))
      } else {
        light = (sphere(x=offset_dist*2,y=offset_dist*2,z=offset_dist*2,
                        radius = widest/2,
                        material = light(intensity=lightintensity))) %>%
          add_object(sphere(x=-offset_dist*2,y=offset_dist*2,z=-offset_dist*2,
                            radius = widest/2,
                            material = light(intensity=lightintensity))) %>%
          add_object(sphere(y=-offset_dist*4,
                            radius=widest/2,
                            material = light(intensity=lightintensity)))
      }
      scene = scene %>%
        add_object(light)

    }
    render_scene(scene = scene,
                 fov = fov, lookfrom = c(0,0,widest*5), ...)
  } else {
    if(is.character(lights)) {
      if(lights != "none") {
        if (lights == "top") {
          light = add_light(directional_light(c(1,1,1)), directional_light(c(1,1,-1)))
        } else {
          light = add_light(directional_light(c(1,1,1)), add_light(directional_light(c(-1,1,-1)),
            directional_light(c(0,-1,0))))
        }
      }
    } else if(is.matrix(lights)) {
      light = lights
    }
    scene = rotate_mesh(scene, angle=angle,order_rotation=order_rotation)
    rasterize_scene(scene = scene, lookat=c(0,0,0),
                    light_info = light,
                    fov = fov, lookfrom = c(0,0,widest*5), ...)
  }
}

#' Render Molecule Model
#'
#' Automatically plots the molecule with a camera position and field of view that includes the full model.
#' For more control over the scene, pass the scene to `rayrender::render_scene()` or `rayvertex::rasterize_scene()` and specify
#' the camera position manually.
#'
#' @param scene Raymesh scene of molecule model.
#' @param width Default `800`. Image width in pixels.
#' @param height Default `800`. Image height in pixels.
#' @param fov Default `NULL`, automatically calculated. Camera field of view.
#' @param lookat Default `NULL`, automatically calculated as `c(0,0,0)`.
#' Camera target point.
#' @param lookfrom Default `NULL`, automatically calculated from the scene
#' bounds. Camera position.
#' @param angle Default `c(0,0,0)`. Degrees to rotate the model around the X, Y, and Z axes. If this
#' is a single number, it is interpreted as `c(0, angle, 0)` and only rotates
#' around the Y axis.
#' @param order_rotation Default `c(1,2,3)`. What order to apply the rotations specified in `angle`.
#' @param lights Default `top`. If `none`, removes all lights. If `bottom`, lights scene with light
#' underneath model. If `both`, adds lights both above and below model. This can also be a matrix of
#' light information generated with `rayvertex`.
#' @param lightintensity Default `80`. Light intensity for pathtraced scenes.
#' @param pathtrace Default `TRUE`. If `TRUE`, convert the raymesh scene to a
#'   `rayrender` raymesh and pathtrace it. Scene-generator rayrender materials
#'   are passed to `rayrender::raymesh_model()` with `override_material = TRUE`
#'   when available. If `FALSE`, rasterize the raymesh scene with `rayvertex`.
#' @param plot Default `TRUE`. If `TRUE`, plot the rendered image to the
#'   current graphics device. If `FALSE`, only return the rendered rayimage
#'   array.
#' @param scene_elements Default `NULL`. Additional backend scene elements to
#'   add after automatic camera calculation and model rotation. For
#'   `pathtrace = TRUE`, elements may be rayrender scenes or rayvertex raymesh
#'   scenes. For `pathtrace = FALSE`, elements must be rayvertex raymesh scenes.
#' @param ... Other arguments to pass to `rayrender::render_scene()` or
#' `rayvertex::rasterize_scene()`. Arguments that only apply to the other
#' backend are ignored. For pathtraced scenes, `background` sets both
#' `backgroundlow` and `backgroundhigh` in `rayrender::render_scene()` and
#' turns on `ambient_light`. For raster scenes, `background` is passed through
#' to `rayvertex::rasterize_scene()`.
#'
#' @return Rendered image
#' @importFrom rayrender render_scene glossy sphere segment add_object light group_objects raymesh_model
#' @importFrom rayvertex rasterize_scene material_list sphere_mesh segment_mesh add_shape directional_light add_light rotate_mesh
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' compact_model = read_pdb(download_pdb("6k16", out_dir = tempdir()))
#' compact_scene = compact_model |>
#'   generate_ribbon_scene(
#'     material = rayrender::diffuse,
#'     show_hetero_atoms = FALSE
#'   )
#'
#' # Start with render_model's default camera, top lighting, and calculated FOV.
#' compact_scene |>
#'   render_model(pathtrace = TRUE, width = 800, height = 800, samples = 128)
#'
#' # Use `plot = FALSE` when you want the image array without drawing it.
#' compact_image = compact_scene |>
#'   render_model(
#'     pathtrace = TRUE,
#'     width = 800,
#'     height = 800,
#'     samples = 128,
#'     plot = FALSE
#'   )
#'
#' # This render uses an explicit camera, scalar Y rotation, a solid
#' # background, weaker top light, and rayrender sampling options. The
#' # camera orientation was extracted interactively in rayrender by 
#' # pressing the "P" key.
#' compact_scene |>
#'   render_model(
#'     pathtrace = TRUE,
#'     width = 800,
#'     height = 800,
#'     fov = 72,
#'     lookat = c(-1.93, 11.99, 8.47),
#'     lookfrom = c(-25.50, 26.52, 14.99),
#'     aperture = 2,
#'     angle = 35,
#'     background = "grey2",
#'     lights = "top",
#'     lightintensity = 50,
#'     samples = 128,
#'     sample_method = "sobol_blue",
#'     clamp_value = 10
#'   )
#'
#' # Vector rotation and a custom rotation order show the same scene from a
#' # different orientation with lights above and below.
#' compact_scene |>
#'   render_model(
#'     pathtrace = TRUE,
#'     angle = c(20, 90, 90),
#'     lights = "both",
#'     fov = 28,
#'     samples = 128,
#'     sample_method = "sobol_blue"
#'   )
#'
#' # Disable render_model's automatic lights and add your own.
#' tmp_exr = tempfile(fileext=".exr")
#' background_light = matrix(0,nrow=500,ncol=1000)
#' background_light[140:160, 250:270] = 700 #key
#' background_light[100:170, 200:300+500] = 10 #fill
#' rayimage::ray_write_image(background_light, tmp_exr)
#' #Generate a key and fill light
#' rayimage::plot_image(background_light)
#' compact_scene |>
#'   render_model(
#'     pathtrace = TRUE,
#'     lights = "none",
#'     background = "grey12",
#'     samples = 128,
#'     fov = 25,
#'     sample_method = "sobol_blue",
#'     rotate_env = 180,
#'     environment_light = tmp_exr
#'   )
#'
#' # Extra rayrender scene elements can be added after the automatic camera and
#' # model rotation are set, which is useful for fixed reference marks.
#' compact_scene |>
#'   render_model(
#'     pathtrace = TRUE,
#'     angle = 35,
#'     scene_elements = rayrender::sphere(
#'       x = -10,
#'       y = -10,
#'       z = -10,
#'       radius = 0.5,
#'       material = rayrender::light(intensity = 15)
#'     ),
#'     samples = 128,
#'     sample_method = "sobol_blue"
#'   )
#'
#' # Raster scenes use rayvertex lighting. This example passes a directional
#' # light matrix through render_model to rasterize_scene().
#' raster_scene = compact_model |>
#'   generate_ribbon_scene()
#' raster_scene |>
#'   render_model(
#'     pathtrace = FALSE,
#'     background = "grey12",
#'     lights = rayvertex::directional_light(c(0.2, 1, -1))
#'   )
render_model = function(
  scene,
  width = 800,
  height = 800,
  fov = NULL,
  lookat = NULL,
  lookfrom = NULL,
  angle = c(0, 0, 0),
  order_rotation = c(1, 2, 3),
  lights = "top",
  lightintensity = 80,
  pathtrace = TRUE,
  plot = TRUE,
  scene_elements = NULL,
  ...
) {
  width = validate_render_dimension(width, "width")
  height = validate_render_dimension(height, "height")
  pathtrace = validate_render_flag(pathtrace, "pathtrace")
  plot = validate_render_flag(plot, "plot")
  render_args = list(...)
  if (length(angle) == 1) {
    angle = c(0, angle, 0)
  }
  mesh_scene = is_raymesh_scene(scene)
  rayrender_scene = is_rayrender_scene(scene)
  if (!mesh_scene && !rayrender_scene) {
    stop("scene must be a raymesh scene")
  }
  if (!pathtrace && !mesh_scene) {
    stop("pathtrace = FALSE requires a raymesh scene")
  }
  if (mesh_scene) {
    mesh_bbox = rayvertex::get_mesh_bbox(scene)
    bbox_x = as.numeric(mesh_bbox[, "x"])
    bbox_y = as.numeric(mesh_bbox[, "y"])
    bbox_z = as.numeric(mesh_bbox[, "z"])
    widest = max(c(abs(bbox_x), abs(bbox_y), abs(bbox_z)))
    max_sphere_radii = 0
    offset_dist = widest + widest / 5
  } else {
    is_not_light = unlist(lapply(scene$material, \(x) x$type)) != "light"

    bbox_x = c()
    bbox_y = c()
    bbox_z = c()
    max_sphere_radii = 0

    scene_model = scene[
      is_not_light &
        (scene$shape == "cylinder" | scene$shape == "sphere"),
    ]
    if (nrow(scene_model) > 0) {
      bbox_x = range(c(bbox_x, scene_model$x), na.rm = TRUE)
      bbox_y = range(c(bbox_y, scene_model$y), na.rm = TRUE)
      bbox_z = range(c(bbox_z, scene_model$z), na.rm = TRUE)
      spheresizes = unlist(lapply(scene$shape_info, \(x) x[[1]]$radius))[
        is_not_light
      ]
      if (length(spheresizes) > 0) {
        max_sphere_radii = max(spheresizes, na.rm = TRUE)
      }
    }

    raymesh_rows = which(is_not_light & scene$shape == "raymesh")
    if (length(raymesh_rows) > 0) {
      for (row_index in raymesh_rows) {
        mesh = scene$shape_info[[row_index]]$mesh_info[[1]]
        mesh_bbox = rayvertex::get_mesh_bbox(mesh)
        mesh_scale = scene$transforms[[row_index]]$scale[[1]]
        if (length(mesh_scale) == 1) {
          mesh_scale = c(mesh_scale, mesh_scale, mesh_scale)
        }
        mesh_x = as.numeric(mesh_bbox[, "x"]) *
          mesh_scale[1] +
          scene$x[row_index]
        mesh_y = as.numeric(mesh_bbox[, "y"]) *
          mesh_scale[2] +
          scene$y[row_index]
        mesh_z = as.numeric(mesh_bbox[, "z"]) *
          mesh_scale[3] +
          scene$z[row_index]
        bbox_x = range(c(bbox_x, mesh_x), na.rm = TRUE)
        bbox_y = range(c(bbox_y, mesh_y), na.rm = TRUE)
        bbox_z = range(c(bbox_z, mesh_z), na.rm = TRUE)
      }
    }

    if (length(bbox_x) == 0) {
      stop("Scene does not contain any renderable objects")
    }

    widest = max(c(abs(bbox_x), abs(bbox_y), abs(bbox_z)))
    offset_dist = widest + widest / 5 + max_sphere_radii
  }
  if (is.null(fov)) {
    fov = atan2(widest + widest / 5 + max_sphere_radii, widest * 5) /
      pi *
      180 *
      2
  }
  if (is.null(lookat)) {
    lookat = c(0, 0, 0)
  } else {
    lookat = validate_camera_vector(lookat, "lookat")
  }
  if (is.null(lookfrom)) {
    lookfrom = c(0, 0, -widest * 5)
  } else {
    lookfrom = validate_camera_vector(lookfrom, "lookfrom")
  }
  if (pathtrace) {
    if (mesh_scene) {
      if (any(angle != 0)) {
        scene = rotate_mesh(
          scene,
          angle = angle,
          order_rotation = order_rotation
        )
      }
      scene = raymesh_scene_to_rayrender_scene(scene)
    } else if (any(angle != 0)) {
      scene = group_objects(
        scene,
        angle = angle,
        order_rotation = order_rotation
      )
    }
    scene = add_render_scene_elements(
      scene = scene,
      scene_elements = scene_elements,
      backend = "pathtrace"
    )
    if (lights != "none") {
      if (lights == "top") {
        light = sphere(
          x = offset_dist * 2,
          y = offset_dist * 2,
          z = -offset_dist * 2,
          radius = widest / 2,
          material = light(intensity = lightintensity, invisible = TRUE)
        ) |>
          add_object(sphere(
            x = -offset_dist * 2,
            y = offset_dist * 2,
            z = offset_dist * 2,
            radius = widest / 2,
            material = light(
              intensity = lightintensity,
              invisible = TRUE
            )
          ))
      } else {
        light = (sphere(
          x = offset_dist * 2,
          y = offset_dist * 2,
          z = -offset_dist * 2,
          radius = widest / 2,
          material = light(intensity = lightintensity, invisible = TRUE)
        )) |>
          add_object(sphere(
            x = -offset_dist * 2,
            y = offset_dist * 2,
            z = offset_dist * 2,
            radius = widest / 2,
            material = light(
              intensity = lightintensity,
              invisible = TRUE
            )
          )) |>
          add_object(sphere(
            y = -offset_dist * 4,
            radius = widest / 2,
            material = light(
              intensity = lightintensity,
              invisible = TRUE
            )
          ))
      }
      scene = scene |>
        add_object(light)
    }
    render_args = translate_pathtrace_render_args(render_args)
    render_args = filter_backend_render_args(render_args, backend = "pathtrace")
    do.call(
      render_scene,
      c(
        list(
          scene = scene,
          width = width,
          height = height,
          fov = fov,
          lookat = lookat,
          lookfrom = lookfrom,
          plot_scene = plot
        ),
        render_args
      )
    )
  } else {
    light = NULL
    if (is.character(lights)) {
      if (lights != "none") {
        if (lights == "top") {
          light = add_light(
            directional_light(c(1, 1, 1)),
            directional_light(c(1, 1, -1))
          )
        } else {
          light = add_light(
            directional_light(c(1, 1, -1)),
            add_light(
              directional_light(c(1, 1, 1)),
              directional_light(c(0, -1, 0))
            )
          )
        }
      }
    } else if (is.matrix(lights)) {
      light = lights
    }
    render_args = filter_backend_render_args(render_args, backend = "raster")
    scene = rotate_mesh(scene, angle = angle, order_rotation = order_rotation)
    scene = add_render_scene_elements(
      scene = scene,
      scene_elements = scene_elements,
      backend = "raster"
    )
    do.call(
      rasterize_scene,
      c(
        list(
          scene = scene,
          width = width,
          height = height,
          lookat = lookat,
          light_info = light,
          fov = fov,
          lookfrom = lookfrom,
          plot = plot
        ),
        render_args
      )
    )
  }
}

#' @keywords internal
validate_render_dimension = function(value, argument_name) {
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(sprintf("%s must be a positive finite number", argument_name))
  }
  return(as.integer(value))
}

#' @keywords internal
validate_camera_vector = function(vector, argument_name) {
  if (length(vector) != 3L || any(!is.finite(vector))) {
    stop(sprintf(
      "%s must be a finite numeric vector of length 3",
      argument_name
    ))
  }
  return(as.numeric(vector))
}

#' @keywords internal
is_raymesh_scene = function(scene) {
  inherits(scene, "ray_mesh") &&
    is.list(scene) &&
    !is.null(scene$vertices) &&
    !is.null(scene$shapes)
}

#' @keywords internal
is_rayrender_scene = function(scene) {
  is.data.frame(scene) &&
    "shape" %in% names(scene) &&
    "material" %in% names(scene) &&
    "shape_info" %in% names(scene)
}

#' @keywords internal
validate_render_flag = function(value, argument_name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be TRUE or FALSE", argument_name))
  }
  return(value)
}

#' @keywords internal
translate_pathtrace_render_args = function(render_args) {
  if ("background" %in% names(render_args)) {
    background = render_args[["background"]]
    render_args[["background"]] = NULL
    render_args[["backgroundlow"]] = background
    render_args[["backgroundhigh"]] = background
    render_args[["ambient_light"]] = TRUE
  }
  render_args[["plot_scene"]] = NULL
  return(render_args)
}

#' @keywords internal
filter_backend_render_args = function(render_args, backend) {
  if (length(render_args) == 0L) {
    return(render_args)
  }

  arg_names = names(render_args)
  if (is.null(arg_names)) {
    return(list())
  }

  allowed_args = switch(
    backend,
    pathtrace = names(formals(rayrender::render_scene)),
    raster = names(formals(rayvertex::rasterize_scene)),
    stop("Unknown render backend")
  )
  keep = nzchar(arg_names) & arg_names %in% allowed_args
  return(render_args[keep])
}

#' @keywords internal
raymesh_scene_to_rayrender_scene = function(scene) {
  shape_groups = split(
    seq_along(scene$shapes),
    vapply(
      seq_along(scene$shapes),
      function(shape_index)
        raymesh_shape_pathtrace_material_key(
          scene,
          shape_index
        ),
      character(1)
    )
  )

  shape_scenes = vector("list", length(shape_groups))
  for (group_index in seq_along(shape_groups)) {
    shape_mesh = subset_raymesh_shapes(scene, shape_groups[[group_index]])
    pathtrace_material = raymesh_shape_pathtrace_material(shape_mesh)
    raymesh_args = list(
      mesh = shape_mesh,
      calculate_consistent_normals = FALSE,
      recalculate_normals = FALSE
    )
    if (!is.null(pathtrace_material)) {
      raymesh_args$override_material = TRUE
      raymesh_args$material = pathtrace_material
    } else {
      raymesh_args$override_material = FALSE
    }
    shape_scenes[[group_index]] = do.call(raymesh_model, raymesh_args)
  }

  Reduce(add_object, shape_scenes)
}

#' @keywords internal
subset_raymesh_shape = function(scene, shape_index) {
  subset_raymesh_shapes(scene, shape_index)
}

#' @keywords internal
subset_raymesh_shapes = function(scene, shape_indices) {
  material_hashes = attr(scene, "material_hashes")
  shape_mesh = list(
    shapes = scene$shapes[shape_indices],
    vertices = scene$vertices[shape_indices],
    texcoords = scene$texcoords[shape_indices],
    normals = scene$normals[shape_indices],
    materials = scene$materials[shape_indices]
  )
  if (
    !is.null(material_hashes) && length(material_hashes) >= max(shape_indices)
  ) {
    attr(shape_mesh, "material_hashes") = material_hashes[shape_indices]
  } else {
    attr(shape_mesh, "material_hashes") = rep("", length(shape_indices))
  }
  class(shape_mesh) = c("ray_mesh", "list")
  return(shape_mesh)
}

#' @keywords internal
raymesh_shape_pathtrace_material = function(shape_mesh) {
  if (length(shape_mesh$materials) == 0L) {
    return(NULL)
  }
  pathtrace_materials = vector("list", length(shape_mesh$materials))
  for (material_index in seq_along(shape_mesh$materials)) {
    if (length(shape_mesh$materials[[material_index]]) != 1L) {
      return(NULL)
    }
    material_ids = shape_mesh$shapes[[material_index]]$material_ids
    if (length(material_ids) == 0L || any(material_ids != 0L)) {
      return(NULL)
    }
    pathtrace_material = attr(
      shape_mesh$materials[[material_index]][[1]],
      "raymolecule_rayrender_material"
    )
    if (!inherits(pathtrace_material, "ray_material")) {
      return(NULL)
    }
    pathtrace_materials[[material_index]] = pathtrace_material
  }

  material_keys = vapply(
    pathtrace_materials,
    function(pathtrace_material) {
      paste(as.character(serialize(pathtrace_material, NULL)), collapse = "")
    },
    character(1)
  )
  if (length(unique(material_keys)) == 1L) {
    return(pathtrace_materials[[1]])
  }
  return(NULL)
}

#' @keywords internal
raymesh_shape_pathtrace_material_key = function(scene, shape_index) {
  pathtrace_material = raymesh_shape_pathtrace_material(
    subset_raymesh_shape(scene, shape_index)
  )
  if (is.null(pathtrace_material)) {
    return("rayvertex")
  }
  paste(as.character(serialize(pathtrace_material, NULL)), collapse = "")
}

#' @keywords internal
add_render_scene_elements = function(scene, scene_elements, backend) {
  elements = normalize_render_scene_elements(scene_elements)
  if (length(elements) == 0L) {
    return(scene)
  }

  for (element in elements) {
    if (identical(backend, "pathtrace")) {
      if (is_raymesh_scene(element)) {
        element = raymesh_scene_to_rayrender_scene(element)
      } else if (!is_rayrender_scene(element)) {
        stop(
          "scene_elements must contain rayrender scenes or rayvertex raymesh scenes"
        )
      }
      scene = add_object(scene, element)
    } else if (identical(backend, "raster")) {
      if (!is_raymesh_scene(element)) {
        stop(
          "scene_elements must contain rayvertex raymesh scenes when pathtrace = FALSE"
        )
      }
      scene = add_shape(scene, element)
    } else {
      stop("Unknown render backend")
    }
  }

  return(scene)
}

#' @keywords internal
normalize_render_scene_elements = function(scene_elements) {
  if (is.null(scene_elements)) {
    return(list())
  }
  if (is_raymesh_scene(scene_elements) || is_rayrender_scene(scene_elements)) {
    return(list(scene_elements))
  }
  if (!is.list(scene_elements) || is.data.frame(scene_elements)) {
    stop(
      "scene_elements must be a scene object or a list of scene objects"
    )
  }

  elements = Filter(Negate(is.null), scene_elements)
  invalid = vapply(
    elements,
    function(element) {
      !is_raymesh_scene(element) && !is_rayrender_scene(element)
    },
    logical(1)
  )
  if (any(invalid)) {
    stop(
      "scene_elements must contain rayrender scenes or rayvertex raymesh scenes"
    )
  }

  return(elements)
}

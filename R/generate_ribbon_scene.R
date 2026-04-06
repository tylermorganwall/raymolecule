#' Build Scene (ribbon)
#'
#' Generates a mesh-based protein ribbon from a PDB model parsed with
#' [read_pdb()]. The ribbon follows a CA-driven centripetal Catmull-Rom spline,
#' uses O atoms to guide orientation, and is emitted as indexed watertight mesh
#' data with UV coordinates. Contiguous `SHEET` annotations are rendered with
#' tapered C-terminal arrowheads.
#'
#' @param model Model extracted from a PDB file with [read_pdb()].
#' @param x Default `0`. X offset, applied after centering.
#' @param y Default `0`. Y offset, applied after centering.
#' @param z Default `0`. Z offset, applied after centering.
#' @param scale Default `1`. Amount to scale the ribbon geometry.
#' @param center Default `TRUE`. Centers the bounding box of the model.
#' @param pathtrace Default `TRUE`. If `FALSE`, a `rayvertex` mesh scene is
#'   returned. If `TRUE`, a `rayrender` scene is returned from the same mesh
#'   data.
#' @param ribbon_width Default `1.6`. Width of the ribbon cross-section.
#' @param ribbon_thickness Default `0.25`. Thickness of the ribbon
#'   cross-section.
#' @param cross_section_resolution Default `24`. Number of perimeter vertices
#'   used to approximate the ribbon cross-section.
#' @param subdivisions Default `8`. Minimum number of spline samples per
#'   residue interval. Longer backbone spans are automatically refined to a
#'   smaller internal step size to avoid visible faceting.
#' @param color_mode Either `"chain"` or `"uv"`. If omitted, single-chain
#'   proteins default to `"uv"` and multi-chain proteins default to `"chain"`.
#' @param chain_colors Optional named vector or named list keyed by chain ID.
#'   Used when `color_mode = "chain"`.
#' @param texture Optional texture path. Used when `color_mode = "uv"`. If
#'   omitted in UV mode, a built-in rainbow gradient is used.
#' @param material Default `rayrender::glossy`. Rayrender material template to
#'   use when `pathtrace = TRUE`.
#' @param material_vertex Default `rayvertex::material_list(type = "phong")`.
#'   Rayvertex material template to use when `pathtrace = FALSE`.
#'
#' @return Rayrender/rayvertex scene containing a ribbon mesh.
#' @export
#'
#' @examples
#' if (file.exists("3nir.pdb")) {
#'   read_pdb("3nir.pdb") |>
#'     generate_ribbon_scene(pathtrace = FALSE) |>
#'     render_model()
#' }
generate_ribbon_scene = function(
	model,
	x = 0,
	y = 0,
	z = 0,
	scale = 1,
	center = TRUE,
	pathtrace = TRUE,
	ribbon_width = 1.6,
	ribbon_thickness = 0.25,
	cross_section_resolution = 24,
	subdivisions = 8,
	color_mode = c("chain", "uv"),
	chain_colors = NULL,
	texture = NULL,
	material = rayrender::glossy,
	material_vertex = rayvertex::material_list(type = "phong")
) {
	if (!identical(model$pdb_type, "pdb") || !is.data.frame(model$residues)) {
		stop("generate_ribbon_scene() requires a PDB model produced by read_pdb()")
	}
	if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
		stop("scale must be a positive number")
	}
	if (missing(color_mode)) {
		protein_chain_ids = unique(model$residues$chain_id[model$residues$has_ca])
		protein_chain_ids = protein_chain_ids[!is.na(protein_chain_ids)]
		color_mode = if (length(protein_chain_ids) <= 1L) {
			"uv"
		} else {
			"chain"
		}
	} else {
		color_mode = match.arg(color_mode)
	}
	if (identical(color_mode, "uv") && is.null(texture)) {
		texture = default_ribbon_texture()
	}
	mesh_data = build_ribbon_mesh(
		residues = model$residues,
		ribbon_width = ribbon_width,
		ribbon_thickness = ribbon_thickness,
		cross_section_resolution = cross_section_resolution,
		subdivisions = subdivisions
	)

	chain_ids = vapply(mesh_data$chains, `[[`, character(1), "chain_id")
	chain_color_lookup = NULL
	if (identical(color_mode, "chain")) {
		chain_color_lookup = resolve_ribbon_chain_colors(chain_ids, chain_colors)
	}

	reference_points = ribbon_reference_points(model)
	center_shift = c(0, 0, 0)
	if (center) {
		center_shift = apply(reference_points, 2, function(values) {
			mean(range(values))
		})
	}
	offset = c(x, y, z)

	mesh_list = vector(mode = "list", length = length(mesh_data$chains))
	for (i in seq_along(mesh_data$chains)) {
		chain_mesh = mesh_data$chains[[i]]
		chain_mesh$vertices = transform_ribbon_vertices(
			vertices = chain_mesh$vertices,
			center_shift = center_shift,
			scale = scale,
			offset = offset
		)

		mesh_material = prepare_ribbon_material(
			pathtrace = pathtrace,
			material = material,
			material_vertex = material_vertex,
			color_mode = color_mode,
			chain_id = chain_mesh$chain_id,
			chain_color = if (is.null(chain_color_lookup)) {
				NULL
			} else {
				chain_color_lookup[chain_mesh$chain_id][[1]]
			},
			texture = texture
		)

		mesh_list[[i]] = rayvertex::construct_mesh(
			vertices = chain_mesh$vertices,
			indices = chain_mesh$indices,
			normals = chain_mesh$normals,
			norm_indices = chain_mesh$norm_indices,
			texcoords = chain_mesh$texcoords,
			tex_indices = chain_mesh$tex_indices,
			material = mesh_material
		)
	}

	mesh_scene = rayvertex::scene_from_list(mesh_list)

	if (!pathtrace) {
		return(mesh_scene)
	}

	return(rayrender::raymesh_model(
		mesh_scene,
		override_material = FALSE,
		calculate_consistent_normals = FALSE,
		recalculate_normals = FALSE
	))
}

#' @keywords internal
prepare_ribbon_material = function(
	pathtrace,
	material,
	material_vertex,
	color_mode,
	chain_id,
	chain_color,
	texture
) {
	if (pathtrace) {
		mesh_material = rayrender_material_to_vertex_material(material)
	} else {
		mesh_material = material_vertex
	}

	if (identical(color_mode, "chain")) {
		mesh_material$diffuse = convert_color(chain_color)
		mesh_material$ambient = convert_color(chain_color)
		mesh_material$ambient_intensity = max(0.2, mesh_material$ambient_intensity)
		mesh_material$diffuse_texname = ""
		mesh_material$ambient_texname = ""
	} else if (!is.null(texture)) {
		mesh_material$diffuse = c(1, 1, 1)
		mesh_material$ambient = c(1, 1, 1)
		mesh_material$ambient_intensity = 1
		mesh_material$diffuse_texname = texture
	}

	return(mesh_material)
}

#' @keywords internal
default_ribbon_texture = local({
	texture_path = NULL

	function() {
		if (!is.null(texture_path) && file.exists(texture_path)) {
			return(texture_path)
		}

		texture_path = file.path(tempdir(), "raymolecule-ribbon-rainbow.png")
		if (!file.exists(texture_path)) {
			palette = grDevices::hcl.colors(256, palette = "Spectral", rev = TRUE)
			grDevices::png(
				filename = texture_path,
				width = length(palette),
				height = 8,
				bg = "transparent"
			)
			on.exit(grDevices::dev.off(), add = TRUE)

			graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
			graphics::plot.new()
			graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
			graphics::rasterImage(
				image = as.raster(matrix(palette, nrow = 1L)),
				xleft = 0,
				ybottom = 0,
				xright = 1,
				ytop = 1,
				interpolate = TRUE
			)
		}

		return(texture_path)
	}
})

#' @keywords internal
rayrender_material_to_vertex_material = function(material) {
	material_info = material()
	material_info = material_info[[1]]
	if (!material_info$type %in% c(7L, 1L, 3L)) {
		stop("material() must be either `glossy`, `diffuse`, or `dielectric`")
	}

	base_color = material_info$properties[[1]]
	mesh_material = rayvertex::material_list(
		diffuse = base_color,
		ambient = base_color,
		ambient_intensity = 0.2,
		type = "phong"
	)

	if (material_info$type == 1L) {
		mesh_material$type = "diffuse"
		mesh_material$specular = c(0, 0, 0)
	} else if (material_info$type == 3L) {
		mesh_material$type = "phong"
		if (length(material_info$properties[[1]]) >= 4L) {
			mesh_material$ior = material_info$properties[[1]][4]
		}
	}

	if (!is.null(material_info$image) && nzchar(material_info$image)) {
		mesh_material$diffuse_texname = material_info$image
	}

	return(mesh_material)
}

#' @keywords internal
resolve_ribbon_chain_colors = function(chain_ids, chain_colors) {
	unique_chain_ids = unique(chain_ids)

	if (is.null(chain_colors)) {
		palette = grDevices::hcl.colors(
			length(unique_chain_ids),
			palette = "Dynamic"
		)
		names(palette) = unique_chain_ids
		return(as.list(palette))
	}

	if (is.atomic(chain_colors) && !is.list(chain_colors)) {
		if (is.null(names(chain_colors)) || any(names(chain_colors) == "")) {
			stop(
				"chain_colors must be a named vector or named list keyed by chain ID"
			)
		}
		chain_colors = as.list(chain_colors)
	}

	if (
		!is.list(chain_colors) ||
			is.null(names(chain_colors)) ||
			any(names(chain_colors) == "")
	) {
		stop("chain_colors must be a named vector or named list keyed by chain ID")
	}

	missing_ids = setdiff(unique_chain_ids, names(chain_colors))
	if (length(missing_ids) > 0L) {
		stop(
			sprintf(
				"Missing chain_colors entries for chain IDs: %s",
				paste(missing_ids, collapse = ", ")
			)
		)
	}

	return(chain_colors[unique_chain_ids])
}

#' @keywords internal
ribbon_reference_points = function(model) {
	if (!is.null(model$atoms) && nrow(model$atoms) > 0L) {
		return(as.matrix(model$atoms[, c("x", "y", "z"), drop = FALSE]))
	}

	ca_rows = model$residues[model$residues$has_ca, , drop = FALSE]
	if (nrow(ca_rows) == 0L) {
		stop("generate_ribbon_scene() requires protein residues with CA atoms")
	}

	return(as.matrix(ca_rows[, c("ca_x", "ca_y", "ca_z"), drop = FALSE]))
}

#' @keywords internal
transform_ribbon_vertices = function(vertices, center_shift, scale, offset) {
	vertices = sweep(vertices, 2, center_shift, FUN = "-")
	vertices = vertices * scale
	vertices = sweep(vertices, 2, offset, FUN = "+")
	return(vertices)
}

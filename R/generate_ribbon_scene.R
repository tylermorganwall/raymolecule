#' Build Scene (ribbon)
#'
#' Generates a mesh-based protein ribbon from a PDB model parsed with
#' [read_pdb()]. Segments with contiguous peptide-plane backbone data are
#' rendered with a peptide-plane cartoon sweep, while incomplete backbone
#' segments fall back to a CA-driven centripetal Catmull-Rom ribbon. The mesh is
#' emitted as indexed watertight geometry with UV coordinates, and eligible
#' `SHEET` annotations are rendered with peptide-plane strand profiles and
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
#' @param show_hetero_atoms Default `TRUE`. If `TRUE`, display non-protein
#'   `HETATM` records as bare spheres alongside the ribbon.
#' @param show_hetero_bonds Default `TRUE`. If `TRUE`, display bonds between
#'   shown hetero atoms as a ball-and-stick overlay alongside the ribbon.
#' @param show_waters Default `FALSE`. If `TRUE`, include water `HETATM`
#'   records when `show_hetero_atoms = TRUE` and `show_hetero_bonds = TRUE`.
#' @param use_vertex_normals Default `FALSE`. If `TRUE`, attach the ribbon
#'   mesh's swept vertex normals. If `FALSE`, raster scenes omit explicit
#'   normals and pathtraced scenes use flat face normals.
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
	material_vertex = rayvertex::material_list(type = "phong"),
	show_hetero_atoms = TRUE,
	show_hetero_bonds = TRUE,
	show_waters = FALSE,
	use_vertex_normals = FALSE
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
	if (!is.logical(use_vertex_normals) || length(use_vertex_normals) != 1L || is.na(use_vertex_normals)) {
		stop("use_vertex_normals must be TRUE or FALSE")
	}
	if (!is.logical(show_hetero_atoms) || length(show_hetero_atoms) != 1L || is.na(show_hetero_atoms)) {
		stop("show_hetero_atoms must be TRUE or FALSE")
	}
	if (!is.logical(show_hetero_bonds) || length(show_hetero_bonds) != 1L || is.na(show_hetero_bonds)) {
		stop("show_hetero_bonds must be TRUE or FALSE")
	}
	if (!is.logical(show_waters) || length(show_waters) != 1L || is.na(show_waters)) {
		stop("show_waters must be TRUE or FALSE")
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
	visible_atoms = select_ribbon_display_atoms(
		model = model,
		show_hetero_atoms = show_hetero_atoms,
		show_waters = show_waters
	)
	visible_bonds = select_ribbon_display_bonds(
		model = model,
		atoms = visible_atoms,
		show_hetero_bonds = show_hetero_bonds
	)

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

		if (use_vertex_normals) {
			mesh_list[[i]] = rayvertex::construct_mesh(
				vertices = chain_mesh$vertices,
				indices = chain_mesh$indices,
				normals = chain_mesh$normals,
				norm_indices = chain_mesh$norm_indices,
				texcoords = chain_mesh$texcoords,
				tex_indices = chain_mesh$tex_indices,
				material = mesh_material
			)
		} else if (pathtrace) {
			flat_normals = build_flat_mesh_normals(
				vertices = chain_mesh$vertices,
				indices = chain_mesh$indices,
				fallback_normals = chain_mesh$normals,
				fallback_norm_indices = chain_mesh$norm_indices
			)
			mesh_list[[i]] = rayvertex::construct_mesh(
				vertices = chain_mesh$vertices,
				indices = chain_mesh$indices,
				normals = flat_normals$normals,
				norm_indices = flat_normals$norm_indices,
				texcoords = chain_mesh$texcoords,
				tex_indices = chain_mesh$tex_indices,
				material = mesh_material
			)
		} else {
			mesh_list[[i]] = rayvertex::construct_mesh(
				vertices = chain_mesh$vertices,
				indices = chain_mesh$indices,
				texcoords = chain_mesh$texcoords,
				tex_indices = chain_mesh$tex_indices,
				material = mesh_material
			)
		}
	}
	if (!pathtrace && nrow(visible_bonds) > 0L) {
		mesh_list = c(
			mesh_list,
			build_ribbon_bond_meshes(
				atoms = visible_atoms,
				bonds = visible_bonds,
				center_shift = center_shift,
				scale = scale,
				offset = offset,
				material_vertex = material_vertex
			)
		)
	}
	if (!pathtrace && nrow(visible_atoms) > 0L) {
		mesh_list = c(
			mesh_list,
			build_ribbon_atom_meshes(
				atoms = visible_atoms,
				center_shift = center_shift,
				scale = scale,
				offset = offset,
				material_vertex = material_vertex
			)
		)
	}

	mesh_scene = rayvertex::scene_from_list(mesh_list)

	if (!pathtrace) {
		return(mesh_scene)
	}

	ribbon_scene = rayrender::raymesh_model(
		mesh_scene,
		override_material = FALSE,
		calculate_consistent_normals = FALSE,
		recalculate_normals = FALSE
	)
	if (nrow(visible_bonds) > 0L) {
		ribbon_scene = rayrender::add_object(
			ribbon_scene,
			build_ribbon_bond_objects(
				atoms = visible_atoms,
				bonds = visible_bonds,
				center_shift = center_shift,
				scale = scale,
				offset = offset,
				material = material
			)
		)
	}
	if (nrow(visible_atoms) > 0L) {
		ribbon_scene = rayrender::add_object(
			ribbon_scene,
			build_ribbon_atom_objects(
				atoms = visible_atoms,
				center_shift = center_shift,
				scale = scale,
				offset = offset,
				material = material
			)
		)
	}

	return(ribbon_scene)
}

#' @keywords internal
build_flat_mesh_normals = function(
	vertices,
	indices,
	fallback_normals = NULL,
	fallback_norm_indices = NULL
) {
	if (!is.matrix(vertices) || ncol(vertices) != 3L) {
		stop("build_flat_mesh_normals() requires an n x 3 vertex matrix")
	}
	if (!is.matrix(indices) || ncol(indices) != 3L) {
		stop("build_flat_mesh_normals() requires an m x 3 index matrix")
	}
	if (is.null(fallback_normals) != is.null(fallback_norm_indices)) {
		stop(
			"build_flat_mesh_normals() requires fallback_normals and fallback_norm_indices together"
		)
	}
	if (!is.null(fallback_normals)) {
		if (!is.matrix(fallback_normals) || ncol(fallback_normals) != 3L) {
			stop("build_flat_mesh_normals() fallback_normals must be an n x 3 matrix")
		}
		if (!is.matrix(fallback_norm_indices) || ncol(fallback_norm_indices) != 3L) {
			stop("build_flat_mesh_normals() fallback_norm_indices must be an m x 3 matrix")
		}
		if (nrow(fallback_norm_indices) != nrow(indices)) {
			stop("build_flat_mesh_normals() fallback_norm_indices must match indices")
		}
	}

	triangle_count = nrow(indices)
	normals = matrix(NA_real_, nrow = triangle_count, ncol = 3L)
	norm_indices = matrix(rep(seq_len(triangle_count) - 1L, each = 3L), ncol = 3L, byrow = TRUE)
	previous_normal = NULL

	for (triangle_index in seq_len(triangle_count)) {
		vertex_ids = indices[triangle_index, ] + 1L
		edge1 = vertices[vertex_ids[2], ] - vertices[vertex_ids[1], ]
		edge2 = vertices[vertex_ids[3], ] - vertices[vertex_ids[1], ]
		normal = cross(edge1, edge2)
		normal_length = vector_length(normal)

		if (is.finite(normal_length) && normal_length > 1e-8) {
			normalized_normal = normal / normal_length
		} else {
			normalized_normal = resolve_degenerate_flat_mesh_normal(
				triangle_index = triangle_index,
				vertices = vertices,
				vertex_ids = vertex_ids,
				fallback_normals = fallback_normals,
				fallback_norm_indices = fallback_norm_indices,
				previous_normal = previous_normal
			)
		}

		normals[triangle_index, ] = normalized_normal
		previous_normal = normalized_normal
	}

	return(list(normals = normals, norm_indices = norm_indices))
}

#' @keywords internal
resolve_degenerate_flat_mesh_normal = function(
	triangle_index,
	vertices,
	vertex_ids,
	fallback_normals = NULL,
	fallback_norm_indices = NULL,
	previous_normal = NULL
) {
	if (!is.null(fallback_normals)) {
		fallback_ids = fallback_norm_indices[triangle_index, ] + 1L
		fallback_normal = colMeans(fallback_normals[fallback_ids, , drop = FALSE])
		fallback_length = vector_length(fallback_normal)
		if (is.finite(fallback_length) && fallback_length > 1e-8) {
			fallback_normal = fallback_normal / fallback_length
			if (!is.null(previous_normal) && sum(fallback_normal * previous_normal) < 0) {
				fallback_normal = -fallback_normal
			}
			return(fallback_normal)
		}
	}

	if (!is.null(previous_normal) && all(is.finite(previous_normal))) {
		previous_length = vector_length(previous_normal)
		if (previous_length > 1e-8) {
			return(previous_normal / previous_length)
		}
	}

	triangle_vertices = vertices[vertex_ids, , drop = FALSE]
	edge_vectors = rbind(
		triangle_vertices[2, ] - triangle_vertices[1, ],
		triangle_vertices[3, ] - triangle_vertices[1, ],
		triangle_vertices[3, ] - triangle_vertices[2, ]
	)
	edge_lengths = sqrt(rowSums(edge_vectors^2))
	if (any(is.finite(edge_lengths) & edge_lengths > 1e-8)) {
		longest_edge = edge_vectors[which.max(edge_lengths), ]
		return(stable_perpendicular(longest_edge))
	}

	return(c(0, 0, 1))
}

#' @keywords internal
select_ribbon_display_atoms = function(
	model,
	show_hetero_atoms,
	show_waters
) {
	if (!show_hetero_atoms || is.null(model$atoms) || nrow(model$atoms) == 0L) {
		return(empty_pdb_atoms())
	}

	atoms = model$atoms[model$atoms$record == "HETATM", , drop = FALSE]
	if (!show_waters) {
		atoms = atoms[!is_water_residue_name(atoms$res_name), , drop = FALSE]
	}
	return(atoms)
}

#' @keywords internal
select_ribbon_display_bonds = function(model, atoms, show_hetero_bonds) {
	if (
		!show_hetero_bonds ||
			nrow(atoms) == 0L ||
			is.null(model$bonds) ||
			nrow(model$bonds) == 0L
	) {
		return(empty_pdb_bonds())
	}

	visible_atom_ids = atoms$index
	bonds = model$bonds[
		model$bonds$from %in% visible_atom_ids &
			model$bonds$to %in% visible_atom_ids,
		,
		drop = FALSE
	]
	if (nrow(bonds) == 0L) {
		return(empty_pdb_bonds())
	}

	bonds = data.frame(
		from = pmin(bonds$from, bonds$to),
		to = pmax(bonds$from, bonds$to),
		number = bonds$number,
		stringsAsFactors = FALSE
	)
	bonds = bonds[bonds$from != bonds$to, , drop = FALSE]
	if (nrow(bonds) == 0L) {
		return(empty_pdb_bonds())
	}

	bonds = stats::aggregate(number ~ from + to, data = bonds, FUN = max)
	bonds = bonds[order(bonds$from, bonds$to), , drop = FALSE]
	rownames(bonds) = NULL
	return(bonds)
}

#' @keywords internal
is_water_residue_name = function(res_name) {
	res_name = toupper(trimws(res_name))
	return(res_name %in% c("HOH", "WAT", "DOD"))
}

#' @keywords internal
normalize_atom_element_symbol = function(element) {
	element = trimws(as.character(element))
	if (!nzchar(element)) {
		return("C")
	}
	if (nchar(element) == 1L) {
		return(toupper(element))
	}
	return(paste0(
		toupper(substr(element, 1L, 1L)),
		tolower(substr(element, 2L, nchar(element)))
	))
}

#' @keywords internal
display_atom_color = function(element) {
	element = normalize_atom_element_symbol(element)
	return(PeriodicTable::atomColor(element))
}

#' @keywords internal
display_atom_radius = function(element) {
	element = normalize_atom_element_symbol(element)
	return((PeriodicTable::mass(element) / 14)^(1 / 3) / 2)
}

#' @keywords internal
ribbon_display_atom_radius = function(element, scale = 0.45) {
	if (
		length(scale) != 1L ||
			!is.finite(scale) ||
			scale <= 0
	) {
		stop("ribbon_display_atom_radius() scale must be a positive number")
	}
	return(display_atom_radius(element) * scale)
}

#' @keywords internal
ribbon_display_bond_radius = function(element1, element2, scale = 0.6) {
	if (
		length(scale) != 1L ||
			!is.finite(scale) ||
			scale <= 0
	) {
		stop("ribbon_display_bond_radius() scale must be a positive number")
	}
	return(min(
		ribbon_display_atom_radius(element1),
		ribbon_display_atom_radius(element2)
	) * scale)
}

#' @keywords internal
transform_scene_points = function(points, center_shift, scale, offset) {
	points = sweep(points, 2, center_shift, FUN = "-")
	points = points * scale
	points = sweep(points, 2, offset, FUN = "+")
	return(points)
}

#' @keywords internal
build_ribbon_bond_meshes = function(
	atoms,
	bonds,
	center_shift,
	scale,
	offset,
	material_vertex
) {
	if (nrow(bonds) == 0L) {
		return(list())
	}

	positions = transform_scene_points(
		points = as.matrix(atoms[, c("x", "y", "z"), drop = FALSE]),
		center_shift = center_shift,
		scale = scale,
		offset = offset
	)
	rownames(positions) = as.character(atoms$index)

	meshes = vector(mode = "list", length = nrow(bonds) * 2L)
	counter = 1L
	for (i in seq_len(nrow(bonds))) {
		from_atom = atoms[match(bonds$from[i], atoms$index), , drop = FALSE]
		to_atom = atoms[match(bonds$to[i], atoms$index), , drop = FALSE]
		start = positions[as.character(bonds$from[i]), ]
		end = positions[as.character(bonds$to[i]), ]
		midpoint = (start + end) / 2
		radius = ribbon_display_bond_radius(from_atom$type, to_atom$type)

		from_material = material_vertex
		from_color = display_atom_color(from_atom$type)
		from_material$diffuse = convert_color(from_color)
		from_material$ambient = convert_color(from_color)
		from_material$ambient_intensity = max(0.3, from_material$ambient_intensity)
		meshes[[counter]] = rayvertex::segment_mesh(
			start = start,
			end = midpoint,
			radius = radius,
			material = from_material
		)
		counter = counter + 1L

		to_material = material_vertex
		to_color = display_atom_color(to_atom$type)
		to_material$diffuse = convert_color(to_color)
		to_material$ambient = convert_color(to_color)
		to_material$ambient_intensity = max(0.3, to_material$ambient_intensity)
		meshes[[counter]] = rayvertex::segment_mesh(
			start = midpoint,
			end = end,
			radius = radius,
			material = to_material
		)
		counter = counter + 1L
	}

	return(meshes)
}

#' @keywords internal
build_ribbon_atom_meshes = function(
	atoms,
	center_shift,
	scale,
	offset,
	material_vertex
) {
	if (nrow(atoms) == 0L) {
		return(list())
	}

	positions = transform_scene_points(
		points = as.matrix(atoms[, c("x", "y", "z"), drop = FALSE]),
		center_shift = center_shift,
		scale = scale,
		offset = offset
	)
	meshes = vector(mode = "list", length = nrow(atoms))
	for (i in seq_len(nrow(atoms))) {
		atom_material = material_vertex
		atom_color = display_atom_color(atoms$type[i])
		atom_material$diffuse = convert_color(atom_color)
		atom_material$ambient = convert_color(atom_color)
		atom_material$ambient_intensity = max(0.3, atom_material$ambient_intensity)
		meshes[[i]] = rayvertex::sphere_mesh(
			position = positions[i, ],
			radius = ribbon_display_atom_radius(atoms$type[i]),
			low_poly = FALSE,
			material = atom_material
		)
	}

	return(meshes)
}

#' @keywords internal
build_ribbon_bond_objects = function(
	atoms,
	bonds,
	center_shift,
	scale,
	offset,
	material
) {
	if (nrow(bonds) == 0L) {
		return(NULL)
	}

	positions = transform_scene_points(
		points = as.matrix(atoms[, c("x", "y", "z"), drop = FALSE]),
		center_shift = center_shift,
		scale = scale,
		offset = offset
	)
	rownames(positions) = as.character(atoms$index)

	objects = vector(mode = "list", length = nrow(bonds) * 2L)
	counter = 1L
	for (i in seq_len(nrow(bonds))) {
		from_atom = atoms[match(bonds$from[i], atoms$index), , drop = FALSE]
		to_atom = atoms[match(bonds$to[i], atoms$index), , drop = FALSE]
		start = positions[as.character(bonds$from[i]), ]
		end = positions[as.character(bonds$to[i]), ]
		midpoint = (start + end) / 2
		radius = ribbon_display_bond_radius(from_atom$type, to_atom$type)

		objects[[counter]] = rayrender::segment(
			start = start,
			end = midpoint,
			radius = radius,
			material = material(color = display_atom_color(from_atom$type))
		)
		counter = counter + 1L
		objects[[counter]] = rayrender::segment(
			start = midpoint,
			end = end,
			radius = radius,
			material = material(color = display_atom_color(to_atom$type))
		)
		counter = counter + 1L
	}

	return(do.call("rbind", objects))
}

#' @keywords internal
build_ribbon_atom_objects = function(
	atoms,
	center_shift,
	scale,
	offset,
	material
) {
	if (nrow(atoms) == 0L) {
		return(NULL)
	}

	positions = transform_scene_points(
		points = as.matrix(atoms[, c("x", "y", "z"), drop = FALSE]),
		center_shift = center_shift,
		scale = scale,
		offset = offset
	)
	objects = vector(mode = "list", length = nrow(atoms))
	for (i in seq_len(nrow(atoms))) {
		objects[[i]] = rayrender::sphere(
			x = positions[i, 1],
			y = positions[i, 2],
			z = positions[i, 3],
			radius = ribbon_display_atom_radius(atoms$type[i]),
			material = material(color = display_atom_color(atoms$type[i]))
		)
	}

	return(do.call("rbind", objects))
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
				image = grDevices::as.raster(matrix(palette, nrow = 1L)),
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

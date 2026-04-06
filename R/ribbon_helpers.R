#' @keywords internal
catmull_rom_chain = function(points) {
	if (!is.matrix(points) || ncol(points) != 3L) {
		stop("catmull_rom_chain() requires an n x 3 numeric matrix")
	}
	if (nrow(points) < 2L) {
		stop("catmull_rom_chain() requires at least two control points")
	}
	if (any(!is.finite(points))) {
		stop("catmull_rom_chain() control points must be finite")
	}

	extended_points = rbind(
		2 * points[1, ] - points[2, ],
		points,
		2 * points[nrow(points), ] - points[nrow(points) - 1L, ]
	)

	segments = vector(mode = "list", length = nrow(points) - 1L)
	for (segment_index in seq_along(segments)) {
		p0 = extended_points[segment_index, ]
		p1 = extended_points[segment_index + 1L, ]
		p2 = extended_points[segment_index + 2L, ]
		p3 = extended_points[segment_index + 3L, ]

		d01 = sqrt(sum((p1 - p0) * (p1 - p0)))
		d12 = sqrt(sum((p2 - p1) * (p2 - p1)))
		d23 = sqrt(sum((p3 - p2) * (p3 - p2)))

		if (any(c(d01, d12, d23) <= 0)) {
			stop(
				"Ribbon generation does not support duplicate consecutive CA positions"
			)
		}

		t0 = 0
		t1 = t0 + d01^0.5
		t2 = t1 + d12^0.5
		t3 = t2 + d23^0.5

		m1 = (t2 - t1) *
			((p1 - p0) / (t1 - t0) - (p2 - p0) / (t2 - t0) + (p2 - p1) / (t2 - t1))

		m2 = (t2 - t1) *
			((p2 - p1) / (t2 - t1) - (p3 - p1) / (t3 - t1) + (p3 - p2) / (t3 - t2))

		segments[[segment_index]] = list(
			p1 = p1,
			p2 = p2,
			m1 = m1,
			m2 = m2
		)
	}

	return(list(points = points, segments = segments))
}

#' @keywords internal
sample_catmull_rom_chain = function(chain, subdivisions) {
	if (length(subdivisions) != 1L || is.na(subdivisions)) {
		stop("subdivisions must be a positive integer")
	}
	subdivisions = as.integer(subdivisions)
	if (subdivisions < 1L) {
		stop("subdivisions must be a positive integer")
	}

	segment_count = length(chain$segments)
	sample_count = segment_count * subdivisions + 1L
	centerline = matrix(NA_real_, nrow = sample_count, ncol = 3L)
	tangent = matrix(NA_real_, nrow = sample_count, ncol = 3L)
	residue_parameter = numeric(sample_count)

	counter = 1L
	for (segment_index in seq_len(segment_count)) {
		s_values = seq(0, 1, length.out = subdivisions + 1L)
		if (segment_index < segment_count) {
			s_values = s_values[-length(s_values)]
		}

		segment = chain$segments[[segment_index]]
		for (s in s_values) {
			centerline[counter, ] = hermite_position(
				p1 = segment$p1,
				p2 = segment$p2,
				m1 = segment$m1,
				m2 = segment$m2,
				s = s
			)

			derivative = hermite_derivative(
				p1 = segment$p1,
				p2 = segment$p2,
				m1 = segment$m1,
				m2 = segment$m2,
				s = s
			)
			tangent[counter, ] = normalize_vector(
				derivative,
				error_message = "Ribbon sampling produced a zero tangent vector"
			)
			residue_parameter[counter] = (segment_index - 1L) + s
			counter = counter + 1L
		}
	}

	return(list(
		centerline = centerline,
		tangent = tangent,
		residue_parameter = residue_parameter
	))
}

#' @keywords internal
build_residue_guides = function(residues, residue_parameter) {
	if (nrow(residues) == 0L) {
		stop("build_residue_guides() requires at least one residue")
	}

	guides = matrix(NA_real_, nrow = length(residue_parameter), ncol = 3L)
	control_guides = matrix(NA_real_, nrow = nrow(residues), ncol = 3L)
	valid_guides = residues$has_ca & residues$has_o

	if (any(valid_guides)) {
		control_guides[valid_guides, ] = cbind(
			residues$o_x[valid_guides] - residues$ca_x[valid_guides],
			residues$o_y[valid_guides] - residues$ca_y[valid_guides],
			residues$o_z[valid_guides] - residues$ca_z[valid_guides]
		)
	}

	tolerance = 1e-10
	for (i in seq_along(residue_parameter)) {
		parameter = residue_parameter[[i]]
		left_index = floor(parameter) + 1L
		right_index = ceiling(parameter) + 1L

		if (abs(parameter - round(parameter)) <= tolerance) {
			residue_index = as.integer(round(parameter)) + 1L
			if (valid_guides[residue_index]) {
				guides[i, ] = control_guides[residue_index, ]
			}
			next
		}

		if (
			left_index >= 1L &&
				right_index <= nrow(residues) &&
				valid_guides[left_index] &&
				valid_guides[right_index]
		) {
			fraction = parameter - floor(parameter)
			guides[i, ] =
				(1 - fraction) *
				control_guides[left_index, ] +
				fraction * control_guides[right_index, ]
		}
	}

	return(guides)
}

#' @keywords internal
build_rotation_minimizing_frame = function(
	centerline,
	tangent,
	guides = NULL,
	chain_id = "",
	residue_parameter = seq_len(nrow(centerline)) - 1
) {
	if (!is.matrix(centerline) || ncol(centerline) != 3L) {
		stop(
			"build_rotation_minimizing_frame() requires an n x 3 centerline matrix"
		)
	}
	if (
		!is.matrix(tangent) ||
			ncol(tangent) != 3L ||
			nrow(tangent) != nrow(centerline)
	) {
		stop(
			"build_rotation_minimizing_frame() requires tangent vectors matching the centerline"
		)
	}
	if (is.null(guides)) {
		guides = matrix(NA_real_, nrow = nrow(centerline), ncol = 3L)
	}
	if (
		!is.matrix(guides) || ncol(guides) != 3L || nrow(guides) != nrow(centerline)
	) {
		stop(
			"build_rotation_minimizing_frame() guide vectors must match the centerline"
		)
	}

	sample_count = nrow(centerline)
	transport_frame = build_parallel_transport_frame(
		centerline = centerline,
		tangent = tangent,
		guides = guides
	)

	twist_angles = compute_guide_twist_angles(
		tangent = tangent,
		guides = guides,
		transport_normal = transport_frame$normal,
		transport_binormal = transport_frame$binormal,
		arc_length = transport_frame$arc_length,
		residue_parameter = residue_parameter
	)

	normal = matrix(NA_real_, nrow = sample_count, ncol = 3L)
	binormal = matrix(NA_real_, nrow = sample_count, ncol = 3L)

	for (i in seq_len(sample_count)) {
		rotated = rotate_frame_about_tangent(
			tangent = tangent[i, ],
			normal = transport_frame$normal[i, ],
			binormal = transport_frame$binormal[i, ],
			angle = twist_angles[i]
		)
		normal[i, ] = rotated$normal
		binormal[i, ] = rotated$binormal
	}

	arc_length = transport_frame$arc_length

	sampled = data.frame(
		x = centerline[, 1],
		y = centerline[, 2],
		z = centerline[, 3],
		tangent_x = tangent[, 1],
		tangent_y = tangent[, 2],
		tangent_z = tangent[, 3],
		normal_x = normal[, 1],
		normal_y = normal[, 2],
		normal_z = normal[, 3],
		binormal_x = binormal[, 1],
		binormal_y = binormal[, 2],
		binormal_z = binormal[, 3],
		arc_length = arc_length,
		chain_id = rep(chain_id, sample_count),
		residue_parameter = residue_parameter,
		stringsAsFactors = FALSE
	)

	return(list(
		centerline = centerline,
		tangent = tangent,
		normal = normal,
		binormal = binormal,
		arc_length = arc_length,
		chain_id = rep(chain_id, sample_count),
		residue_parameter = residue_parameter,
		sampled = sampled
	))
}

#' @keywords internal
build_parallel_transport_frame = function(centerline, tangent, guides = NULL) {
	sample_count = nrow(tangent)
	normal = matrix(NA_real_, nrow = sample_count, ncol = 3L)
	binormal = matrix(NA_real_, nrow = sample_count, ncol = 3L)

	if (is.null(guides)) {
		guides = matrix(NA_real_, nrow = sample_count, ncol = 3L)
	}

	seed_normal = NULL
	for (i in seq_len(sample_count)) {
		if (all(is.finite(guides[i, ]))) {
			projected = project_to_plane(guides[i, ], tangent[i, ])
			if (vector_length(projected) > 1e-8) {
				seed_normal = normalize_vector(
					projected,
					error_message = "Failed to initialize the ribbon frame from guide atoms"
				)
				break
			}
		}
	}

	if (is.null(seed_normal)) {
		seed_normal = stable_perpendicular(tangent[1, ])
	}

	normal[1, ] = normalize_vector(
		project_to_plane(seed_normal, tangent[1, ]),
		error_message = "Failed to initialize the ribbon frame"
	)
	binormal[1, ] = normalize_vector(
		cross(tangent[1, ], normal[1, ]),
		error_message = "Failed to initialize the ribbon binormal"
	)
	normal[1, ] = normalize_vector(
		cross(binormal[1, ], tangent[1, ]),
		error_message = "Failed to orthogonalize the ribbon frame"
	)

	for (i in 2:sample_count) {
		normal[i, ] = parallel_transport_normal(
			previous_normal = normal[i - 1L, ],
			previous_tangent = tangent[i - 1L, ],
			current_tangent = tangent[i, ]
		)
		binormal[i, ] = normalize_vector(
			cross(tangent[i, ], normal[i, ]),
			error_message = "Failed to compute a ribbon binormal"
		)
		normal[i, ] = normalize_vector(
			cross(binormal[i, ], tangent[i, ]),
			error_message = "Failed to compute a ribbon normal"
		)
	}

	arc_length = numeric(sample_count)
	if (sample_count > 1L) {
		segment_lengths = sqrt(rowSums(
			(centerline[-1, , drop = FALSE] -
				centerline[-sample_count, , drop = FALSE])^2
		))
		arc_length[-1] = cumsum(segment_lengths)
	}
	return(list(
		normal = normal,
		binormal = binormal,
		arc_length = arc_length
	))
}

#' @keywords internal
compute_guide_twist_angles = function(
	tangent,
	guides,
	transport_normal,
	transport_binormal,
	arc_length,
	residue_parameter
) {
	sample_count = nrow(tangent)
	tolerance = 1e-8
	knot_mask = abs(residue_parameter - round(residue_parameter)) <= tolerance
	knot_indices = which(knot_mask)

	if (length(knot_indices) == 0L) {
		return(rep(0, sample_count))
	}

	knot_angles = rep(NA_real_, length(knot_indices))
	previous_angle = NULL

	for (k in seq_along(knot_indices)) {
		index = knot_indices[[k]]
		guide = guides[index, ]
		if (!all(is.finite(guide))) {
			next
		}

		projected = project_to_plane(guide, tangent[index, ])
		if (vector_length(projected) <= 1e-8) {
			next
		}
		projected = normalize_vector(
			projected,
			error_message = "Failed to orient the ribbon frame from guide atoms"
		)

		raw_angle = atan2(
			sum(projected * transport_binormal[index, ]),
			sum(projected * transport_normal[index, ])
		)
		knot_angles[[k]] = choose_continuous_axial_twist_angle(
			angle = raw_angle,
			previous_angle = previous_angle
		)
		previous_angle = knot_angles[[k]]
	}

	valid_knots = !is.na(knot_angles)
	if (!any(valid_knots)) {
		return(rep(0, sample_count))
	}

  knot_positions = arc_length[knot_indices[valid_knots]]
  return(interpolate_twist_angles(
    positions = arc_length,
    knot_positions = knot_positions,
    knot_angles = knot_angles[valid_knots]
  ))
}

#' @keywords internal
choose_continuous_twist_angle = function(angle, previous_angle = NULL) {
  if (is.null(previous_angle)) {
    return(angle)
  }
  turn_count = round((previous_angle - angle) / (2 * pi))
  return(angle + 2 * pi * turn_count)
}

#' @keywords internal
choose_continuous_axial_twist_angle = function(angle, previous_angle = NULL) {
	directed_angle = choose_continuous_twist_angle(
		angle = angle,
		previous_angle = previous_angle
	)
	flipped_angle = choose_continuous_twist_angle(
		angle = angle + pi,
		previous_angle = previous_angle
	)

	if (is.null(previous_angle)) {
		if (abs(flipped_angle) < abs(directed_angle)) {
			return(flipped_angle)
		}
		return(directed_angle)
	}

	if (abs(flipped_angle - previous_angle) < abs(directed_angle - previous_angle)) {
		return(flipped_angle)
	}

	return(directed_angle)
}

#' @keywords internal
interpolate_twist_angles = function(positions, knot_positions, knot_angles) {
  if (
    length(knot_positions) != length(knot_angles) ||
      length(knot_positions) == 0L
  ) {
    stop("interpolate_twist_angles() requires matching non-empty knot positions and angles")
  }
  if (any(!is.finite(positions)) || any(!is.finite(knot_positions)) || any(!is.finite(knot_angles))) {
    stop("interpolate_twist_angles() requires finite positions and angles")
  }
  if (any(diff(knot_positions) <= 0)) {
    stop("interpolate_twist_angles() requires strictly increasing knot positions")
  }

  if (length(knot_positions) == 1L) {
    return(rep(knot_angles, length(positions)))
  }

  if (length(knot_positions) == 2L) {
    return(stats::approx(
      x = knot_positions,
      y = knot_angles,
      xout = positions,
      method = "linear",
      rule = 2
    )$y)
  }

  spline = stats::splinefun(
    x = knot_positions,
    y = knot_angles,
    method = "natural"
  )
  clamped_positions = pmin(
    pmax(positions, knot_positions[1]),
    knot_positions[length(knot_positions)]
  )
  return(as.numeric(spline(clamped_positions)))
}

#' @keywords internal
smooth_twist_angles = function(positions, angles, lambda = 2.5) {
	if (length(positions) != length(angles)) {
		stop(
			"smooth_twist_angles() positions and angles must have matching lengths"
		)
	}
	if (length(angles) < 3L) {
		return(angles)
	}
	if (any(!is.finite(positions)) || any(!is.finite(angles))) {
		stop("smooth_twist_angles() requires finite positions and angles")
	}
	if (any(diff(positions) <= 0)) {
		stop("smooth_twist_angles() requires strictly increasing positions")
	}

	count = length(angles)
	second_difference = matrix(0, nrow = count - 2L, ncol = count)
	for (i in seq_len(nrow(second_difference))) {
		second_difference[i, i:(i + 2L)] = c(1, -2, 1)
	}

	data_weight = diag(count)
	data_weight[1, 1] = 2
	data_weight[count, count] = 2
	system_matrix = data_weight + lambda * crossprod(second_difference)
	right_hand_side = data_weight %*% angles

	return(as.numeric(solve(system_matrix, right_hand_side)))
}

#' @keywords internal
rotate_frame_about_tangent = function(tangent, normal, binormal, angle) {
	rotated_normal = cos(angle) * normal + sin(angle) * binormal
	rotated_binormal = normalize_vector(
		cross(tangent, rotated_normal),
		error_message = "Failed to rotate the ribbon frame"
	)
	rotated_normal = normalize_vector(
		cross(rotated_binormal, tangent),
		error_message = "Failed to rotate the ribbon frame"
	)
	return(list(normal = rotated_normal, binormal = rotated_binormal))
}

#' @keywords internal
build_ribbon_dimensions = function(
	residues,
	residue_parameter,
	ribbon_width,
	ribbon_thickness,
	arrow_length = 2,
	arrow_width_scale = 1.75,
	arrow_tip_scale = 0.15,
	arrow_recovery_length = 0.75,
	arrow_base_fraction = 0.65
) {
	if (!is.data.frame(residues)) {
		stop("build_ribbon_dimensions() requires a residue data frame")
	}
	if (length(residue_parameter) == 0L) {
		stop("build_ribbon_dimensions() requires sampled residue parameters")
	}

	widths = rep(ribbon_width, length(residue_parameter))
	thicknesses = rep(ribbon_thickness, length(residue_parameter))
	sheet_runs = find_sheet_runs(residues)

	if (length(sheet_runs) == 0L) {
		return(list(width = widths, thickness = thicknesses))
	}

		max_parameter = max(residue_parameter)
	for (run in sheet_runs) {
		tip_parameter = run$end - 1
		base_parameter = max(run$start - 1, tip_parameter - arrow_length)
		head_extent = tip_parameter - base_parameter
		if (head_extent <= 1e-8) {
			tip_mask = abs(residue_parameter - tip_parameter) <= 1e-8
			widths[tip_mask] = ribbon_width * arrow_tip_scale
		} else {
			head_mask =
				residue_parameter >= base_parameter &
				residue_parameter <= tip_parameter
			u = (residue_parameter[head_mask] - base_parameter) / head_extent
			head_scale = numeric(length(u))

			leading = u <= arrow_base_fraction
			if (any(leading)) {
				head_scale[leading] = 1 + (arrow_width_scale - 1) *
					smoothstep01(u[leading] / arrow_base_fraction)
			}
			if (any(!leading)) {
				trailing_u = (u[!leading] - arrow_base_fraction) / (1 - arrow_base_fraction)
				head_scale[!leading] = arrow_width_scale +
					(arrow_tip_scale - arrow_width_scale) * smoothstep01(trailing_u)
			}
				widths[head_mask] = ribbon_width * head_scale
			}

			recovery_end = min(max_parameter, tip_parameter + arrow_recovery_length)
			if (recovery_end > tip_parameter) {
			recovery_mask =
				residue_parameter > tip_parameter &
				residue_parameter <= recovery_end
			recovery_u =
				(residue_parameter[recovery_mask] - tip_parameter) /
				(recovery_end - tip_parameter)
			recovery_scale = arrow_tip_scale +
				(1 - arrow_tip_scale) * smoothstep01(recovery_u)
				widths[recovery_mask] = ribbon_width * recovery_scale
			}
		}

	return(list(width = widths, thickness = thicknesses))
}

#' @keywords internal
find_sheet_runs = function(residues) {
	if (nrow(residues) == 0L) {
		return(list())
	}

	sheet_mask = residues$ss_class == "sheet"
	sheet_keys = ifelse(
		is.na(residues$sheet_id),
		"",
		paste(
			residues$sheet_id,
			if ("sheet_strand" %in% names(residues)) residues$sheet_strand else NA_integer_,
			sep = ":"
		)
	)
	runs = list()
	run_start = NA_integer_
	run_key = NA_character_

	for (i in seq_len(nrow(residues))) {
		if (!sheet_mask[i]) {
			if (!is.na(run_start)) {
				runs[[length(runs) + 1L]] = list(start = run_start, end = i - 1L, key = run_key)
				run_start = NA_integer_
				run_key = NA_character_
			}
			next
		}

		if (is.na(run_start)) {
			run_start = i
			run_key = sheet_keys[i]
			next
		}

		if (!identical(sheet_keys[i], run_key)) {
			runs[[length(runs) + 1L]] = list(start = run_start, end = i - 1L, key = run_key)
			run_start = i
			run_key = sheet_keys[i]
		}
	}

	if (!is.na(run_start)) {
		runs[[length(runs) + 1L]] = list(
			start = run_start,
			end = nrow(residues),
			key = run_key
		)
	}

	return(runs)
}

#' @keywords internal
smoothstep01 = function(x) {
	clamped = pmax(0, pmin(1, x))
	return(clamped * clamped * (3 - 2 * clamped))
}

#' @keywords internal
build_ribbon_mesh = function(
	residues,
	ribbon_width = 1.6,
	ribbon_thickness = 0.25,
	cross_section_resolution = 24,
	subdivisions = 8
) {
	if (!is.data.frame(residues)) {
		stop("build_ribbon_mesh() requires a residue data frame")
	}
	if (nrow(residues) == 0L) {
		stop("Ribbon generation requires protein residues with CA atoms")
	}
	if (
		length(ribbon_width) != 1L || !is.finite(ribbon_width) || ribbon_width <= 0
	) {
		stop("ribbon_width must be a positive number")
	}
	if (
		length(ribbon_thickness) != 1L ||
			!is.finite(ribbon_thickness) ||
			ribbon_thickness <= 0
	) {
		stop("ribbon_thickness must be a positive number")
	}
	if (length(subdivisions) != 1L || is.na(subdivisions)) {
		stop("subdivisions must be a positive integer")
	}
	subdivisions = as.integer(subdivisions)
	if (subdivisions < 1L) {
		stop("subdivisions must be a positive integer")
	}
	cross_section_resolution = normalize_cross_section_resolution(
		cross_section_resolution
	)

	segment_rows = split_ribbon_segments(residues)
	if (length(segment_rows) == 0L) {
		stop("Ribbon generation requires at least two CA atoms in a protein chain")
	}

	chain_meshes = vector(mode = "list", length = length(segment_rows))
	sampled_paths = vector(mode = "list", length = length(segment_rows))

	for (segment_index in seq_along(segment_rows)) {
		chain_residues = residues[segment_rows[[segment_index]], , drop = FALSE]
		control_points = cbind(
			chain_residues$ca_x,
			chain_residues$ca_y,
			chain_residues$ca_z
		)
		chain = catmull_rom_chain(control_points)
		effective_subdivisions = refine_ribbon_subdivisions(
			control_points = control_points,
			subdivisions = subdivisions
		)
		sampled_chain = sample_catmull_rom_chain(
			chain = chain,
			subdivisions = effective_subdivisions
		)
		guides = build_residue_guides(
			residues = chain_residues,
			residue_parameter = sampled_chain$residue_parameter
		)

		frame = build_rotation_minimizing_frame(
			centerline = sampled_chain$centerline,
			tangent = sampled_chain$tangent,
			guides = guides,
			chain_id = chain_residues$chain_id[1],
			residue_parameter = sampled_chain$residue_parameter
		)
		ribbon_dimensions = build_ribbon_dimensions(
			residues = chain_residues,
			residue_parameter = frame$residue_parameter,
			ribbon_width = ribbon_width,
			ribbon_thickness = ribbon_thickness
		)
		frame$sampled$ribbon_width = ribbon_dimensions$width
		frame$sampled$ribbon_thickness = ribbon_dimensions$thickness

		chain_meshes[[segment_index]] = build_single_ribbon_chain_mesh(
			frame = frame,
			ribbon_widths = ribbon_dimensions$width,
			ribbon_thicknesses = ribbon_dimensions$thickness,
			cross_section_resolution = cross_section_resolution,
			chain_id = chain_residues$chain_id[1],
			chain_index = segment_index
		)
		sampled_paths[[segment_index]] = frame
	}

	return(list(chains = chain_meshes, sampled_paths = sampled_paths))
}

#' @keywords internal
refine_ribbon_subdivisions = function(
	control_points,
	subdivisions,
	target_step = 0.15
) {
	if (!is.matrix(control_points) || ncol(control_points) != 3L) {
		stop("refine_ribbon_subdivisions() requires an n x 3 control point matrix")
	}
	if (nrow(control_points) < 2L) {
		stop("refine_ribbon_subdivisions() requires at least two control points")
	}

	segment_lengths = sqrt(rowSums(diff(control_points)^2))
	required_subdivisions = ceiling(max(segment_lengths) / target_step)
	return(max(subdivisions, required_subdivisions))
}

#' @keywords internal
build_ribbon_ring = function(
	center,
	normal,
	binormal,
	profile
) {
	vertices =
		matrix(
			rep(center, each = nrow(profile$vertices)),
			ncol = 3L,
			byrow = FALSE
		) +
		profile$vertices[, 1] *
			matrix(rep(normal, each = nrow(profile$vertices)), ncol = 3L) +
		profile$vertices[, 2] *
			matrix(rep(binormal, each = nrow(profile$vertices)), ncol = 3L)

	normals =
		profile$normals[, 1] *
		matrix(rep(normal, each = nrow(profile$normals)), ncol = 3L) +
		profile$normals[, 2] *
			matrix(rep(binormal, each = nrow(profile$normals)), ncol = 3L)

	return(list(
		vertices = vertices,
		normals = normals,
		v = profile$v
	))
}

#' @keywords internal
compute_ribbon_surface_normals = function(
	vertices,
	centerline,
	ring_count,
	ring_size,
	fallback_normals = NULL
) {
	if (!is.matrix(vertices) || ncol(vertices) != 3L) {
		stop("compute_ribbon_surface_normals() requires an n x 3 vertex matrix")
	}
	if (!is.matrix(centerline) || ncol(centerline) != 3L || nrow(centerline) != ring_count) {
		stop("compute_ribbon_surface_normals() requires centerline rows matching ring_count")
	}
	if (nrow(vertices) != ring_count * ring_size) {
		stop("compute_ribbon_surface_normals() requires vertices matching ring_count * ring_size")
	}
	if (!is.null(fallback_normals)) {
		if (
			!is.matrix(fallback_normals) ||
				ncol(fallback_normals) != 3L ||
				nrow(fallback_normals) != nrow(vertices)
		) {
			stop("compute_ribbon_surface_normals() fallback_normals must match vertices")
		}
	}

	normals = matrix(NA_real_, nrow = nrow(vertices), ncol = 3L)

	for (ring_index in seq_len(ring_count)) {
		prev_ring = max(1L, ring_index - 1L)
		next_ring = min(ring_count, ring_index + 1L)

		for (vertex_index in seq_len(ring_size)) {
			current_row = (ring_index - 1L) * ring_size + vertex_index
			prev_vertex_index = if (vertex_index == 1L) ring_size else vertex_index - 1L
			next_vertex_index = if (vertex_index == ring_size) 1L else vertex_index + 1L
			prev_row = (ring_index - 1L) * ring_size + prev_vertex_index
			next_row = (ring_index - 1L) * ring_size + next_vertex_index
			prev_ring_row = (prev_ring - 1L) * ring_size + vertex_index
			next_ring_row = (next_ring - 1L) * ring_size + vertex_index

			around = vertices[next_row, ] - vertices[prev_row, ]
			along = vertices[next_ring_row, ] - vertices[prev_ring_row, ]
			normal = cross(around, along)

			if (vector_length(normal) <= 1e-8) {
				normal = cross(along, around)
			}
			if (vector_length(normal) <= 1e-8 && !is.null(fallback_normals)) {
				normal = fallback_normals[current_row, ]
			}
			if (vector_length(normal) <= 1e-8) {
				stop("Failed to compute swept ribbon surface normals")
			}

			normal = normalize_vector(
				normal,
				error_message = "Failed to compute swept ribbon surface normals"
			)
			reference = vertices[current_row, ] - centerline[ring_index, ]
			if (sum(normal * reference) < 0) {
				normal = -normal
			}
			normals[current_row, ] = normal
		}
	}

	return(normals)
}

connect_ribbon_rings = function(vertices, ring_count, ring_size) {
	triangle_count = (ring_count - 1L) * ring_size * 2L
	indices = matrix(0L, nrow = triangle_count, ncol = 3L)
	norm_indices = matrix(0L, nrow = triangle_count, ncol = 3L)
	tex_indices = matrix(0L, nrow = triangle_count, ncol = 3L)

	counter = 1L
	for (ring_index in 0:(ring_count - 2L)) {
		vertex_a = ring_index * ring_size
		vertex_b = (ring_index + 1L) * ring_size
		tex_a = ring_index * (ring_size + 1L)
		tex_b = (ring_index + 1L) * (ring_size + 1L)

		for (vertex_index in 0:(ring_size - 1L)) {
			next_index = (vertex_index + 1L) %% ring_size
			seam = vertex_index == ring_size - 1L

			a = vertex_a + vertex_index
			b = vertex_a + next_index
			c = vertex_b + next_index
			d = vertex_b + vertex_index

			diagonal_ac = sum((vertices[c + 1L, ] - vertices[a + 1L, ])^2)
			diagonal_bd = sum((vertices[d + 1L, ] - vertices[b + 1L, ])^2)

			b_tex = if (seam) tex_a + ring_size else tex_a + next_index
			c_tex = if (seam) tex_b + ring_size else tex_b + next_index

			if (diagonal_ac <= diagonal_bd) {
				indices[counter, ] = c(a, b, c)
				indices[counter + 1L, ] = c(a, c, d)
				norm_indices[counter, ] = c(a, b, c)
				norm_indices[counter + 1L, ] = c(a, c, d)
				tex_indices[counter, ] = c(tex_a + vertex_index, b_tex, c_tex)
				tex_indices[counter + 1L, ] = c(
					tex_a + vertex_index,
					c_tex,
					tex_b + vertex_index
				)
			} else {
				indices[counter, ] = c(a, b, d)
				indices[counter + 1L, ] = c(b, c, d)
				norm_indices[counter, ] = c(a, b, d)
				norm_indices[counter + 1L, ] = c(b, c, d)
				tex_indices[counter, ] = c(
					tex_a + vertex_index,
					b_tex,
					tex_b + vertex_index
				)
				tex_indices[counter + 1L, ] = c(b_tex, c_tex, tex_b + vertex_index)
			}
			counter = counter + 2L
		}
	}

	return(list(
		indices = indices,
		norm_indices = norm_indices,
		tex_indices = tex_indices
	))
}

#' @keywords internal
cap_ribbon_ends = function(ring_count, ring_size) {
	triangle_count = ring_size * 2L
	indices = matrix(0L, nrow = triangle_count, ncol = 3L)
	norm_indices = matrix(0L, nrow = triangle_count, ncol = 3L)
	tex_indices = matrix(0L, nrow = triangle_count, ncol = 3L)

	start_center = ring_count * ring_size
	end_center = start_center + 1L
	start_normal = ring_count * ring_size
	end_normal = start_normal + 1L
	start_tex = ring_count * (ring_size + 1L)
	end_tex = start_tex + 1L
	last_ring_vertex = (ring_count - 1L) * ring_size
	last_ring_tex = (ring_count - 1L) * (ring_size + 1L)

	counter = 1L
	for (vertex_index in 0:(ring_size - 1L)) {
		next_index = (vertex_index + 1L) %% ring_size
		start_next_tex = if (vertex_index == ring_size - 1L) {
			ring_size
		} else {
			next_index
		}
		end_next_tex = if (vertex_index == ring_size - 1L) {
			last_ring_tex + ring_size
		} else {
			last_ring_tex + next_index
		}

		indices[counter, ] = c(start_center, next_index, vertex_index)
		norm_indices[counter, ] = c(start_normal, start_normal, start_normal)
		tex_indices[counter, ] = c(start_tex, start_next_tex, vertex_index)
		counter = counter + 1L

		indices[counter, ] = c(
			end_center,
			last_ring_vertex + vertex_index,
			last_ring_vertex + next_index
		)
		norm_indices[counter, ] = c(end_normal, end_normal, end_normal)
		tex_indices[counter, ] = c(
			end_tex,
			last_ring_tex + vertex_index,
			end_next_tex
		)
		counter = counter + 1L
	}

	return(list(
		indices = indices,
		norm_indices = norm_indices,
		tex_indices = tex_indices
	))
}

#' @keywords internal
build_single_ribbon_chain_mesh = function(
	frame,
	ribbon_widths,
	ribbon_thicknesses,
	cross_section_resolution,
	chain_id,
	chain_index
) {
	profile = ribbon_cross_section_profile(
		ribbon_width = ribbon_widths[1],
		ribbon_thickness = ribbon_thicknesses[1],
		cross_section_resolution = cross_section_resolution
	)
	ring_count = nrow(frame$centerline)
	ring = build_ribbon_ring(
		center = frame$centerline[1, ],
		normal = frame$normal[1, ],
		binormal = frame$binormal[1, ],
		profile = profile
	)
	ring_size = nrow(ring$vertices)

	vertex_count = ring_count * ring_size + 2L
	normal_count = ring_count * ring_size + 2L
	texcoord_count = ring_count * (ring_size + 1L) + 2L

	vertices = matrix(NA_real_, nrow = vertex_count, ncol = 3L)
	normals = matrix(NA_real_, nrow = normal_count, ncol = 3L)
	texcoords = matrix(NA_real_, nrow = texcoord_count, ncol = 2L)

	total_length = frame$arc_length[ring_count]
	if (total_length <= 0) {
		u_values = rep(0, ring_count)
	} else {
		u_values = frame$arc_length / total_length
	}

	for (ring_index in seq_len(ring_count)) {
		profile = ribbon_cross_section_profile(
			ribbon_width = ribbon_widths[ring_index],
			ribbon_thickness = ribbon_thicknesses[ring_index],
			cross_section_resolution = cross_section_resolution
		)
		ring = build_ribbon_ring(
			center = frame$centerline[ring_index, ],
			normal = frame$normal[ring_index, ],
			binormal = frame$binormal[ring_index, ],
			profile = profile
		)

		vertex_rows = ((ring_index - 1L) * ring_size + 1L):(ring_index * ring_size)
		tex_rows = ((ring_index - 1L) * (ring_size + 1L) + 1L):(ring_index *
			(ring_size + 1L))

		vertices[vertex_rows, ] = ring$vertices
		normals[vertex_rows, ] = ring$normals
		texcoords[tex_rows[-length(tex_rows)], ] = cbind(
			rep(u_values[ring_index], ring_size),
			ring$v
		)
		texcoords[tex_rows[length(tex_rows)], ] = c(u_values[ring_index], 1)
	}

	side_vertex_count = ring_count * ring_size
	normals[seq_len(side_vertex_count), ] = compute_ribbon_surface_normals(
		vertices = vertices[seq_len(side_vertex_count), , drop = FALSE],
		centerline = frame$centerline,
		ring_count = ring_count,
		ring_size = ring_size,
		fallback_normals = normals[seq_len(side_vertex_count), , drop = FALSE]
	)

	start_center_index = ring_count * ring_size + 1L
	end_center_index = start_center_index + 1L
	vertices[start_center_index, ] = frame$centerline[1, ]
	vertices[end_center_index, ] = frame$centerline[ring_count, ]
	normals[start_center_index, ] = -frame$tangent[1, ]
	normals[end_center_index, ] = frame$tangent[ring_count, ]

	start_center_tex = ring_count * (ring_size + 1L) + 1L
	end_center_tex = start_center_tex + 1L
	texcoords[start_center_tex, ] = c(0, 0.5)
	texcoords[end_center_tex, ] = c(1, 0.5)

	side_data = connect_ribbon_rings(
		vertices = vertices[seq_len(ring_count * ring_size), , drop = FALSE],
		ring_count = ring_count,
		ring_size = ring_size
	)
	cap_data = cap_ribbon_ends(ring_count = ring_count, ring_size = ring_size)

	indices = rbind(side_data$indices, cap_data$indices)
	norm_indices = rbind(side_data$norm_indices, cap_data$norm_indices)
	tex_indices = rbind(side_data$tex_indices, cap_data$tex_indices)

	return(list(
		chain_id = chain_id,
		chain_index = chain_index,
		vertices = vertices,
		indices = indices,
		normals = normals,
		norm_indices = norm_indices,
		texcoords = texcoords,
		tex_indices = tex_indices,
		sampled = frame
	))
}

#' @keywords internal
split_ribbon_segments = function(residues) {
	if (nrow(residues) == 0L) {
		return(list())
	}

	starts = c()
	ends = c()
	run_start = NA_integer_

	for (i in seq_len(nrow(residues))) {
		if (!residues$has_ca[i]) {
			if (!is.na(run_start)) {
				starts = c(starts, run_start)
				ends = c(ends, i - 1L)
				run_start = NA_integer_
			}
			next
		}

		if (
			is.na(run_start) ||
				residues$chain_break_before[i] ||
				(i > 1L && !residues$has_ca[i - 1L])
		) {
			if (!is.na(run_start)) {
				starts = c(starts, run_start)
				ends = c(ends, i - 1L)
			}
			run_start = i
		}

		if (
			residues$chain_break_after[i] ||
				i == nrow(residues) ||
				!residues$has_ca[min(i + 1L, nrow(residues))]
		) {
			starts = c(starts, run_start)
			ends = c(ends, i)
			run_start = NA_integer_
		}
	}

	segments = list()
	segment_counter = 1L
	for (i in seq_along(starts)) {
		if ((ends[[i]] - starts[[i]] + 1L) < 2L) {
			next
		}
		segments[[segment_counter]] = starts[[i]]:ends[[i]]
		segment_counter = segment_counter + 1L
	}

	return(segments)
}

#' @keywords internal
ribbon_cross_section_profile = function(
	ribbon_width,
	ribbon_thickness,
	cross_section_resolution = 24
) {
	cross_section_resolution = normalize_cross_section_resolution(
		cross_section_resolution
	)
	half_width = ribbon_width / 2
	half_thickness = ribbon_thickness / 2
	shape_exponent = 4
	theta = seq(0, 2 * pi, length.out = cross_section_resolution + 1L)
	theta = theta[-length(theta)] - (3 * pi / 4)

	vertices = cbind(
		half_width * sign(cos(theta)) * abs(cos(theta))^(2 / shape_exponent),
		half_thickness * sign(sin(theta)) * abs(sin(theta))^(2 / shape_exponent)
	)

	if (polygon_area(vertices) < 0) {
		vertices = vertices[nrow(vertices):1, , drop = FALSE]
	}

	edge_vectors = vertices[c(2:nrow(vertices), 1L), , drop = FALSE] - vertices
	edge_lengths = sqrt(rowSums(edge_vectors^2))
	perimeter = sum(edge_lengths)
	cumulative = c(0, cumsum(edge_lengths))
	v = cumulative[seq_len(nrow(vertices))] / perimeter

	normals = matrix(NA_real_, nrow = nrow(vertices), ncol = 2L)
	previous_index = c(nrow(vertices), seq_len(nrow(vertices) - 1L))
	next_index = c(2:nrow(vertices), 1L)

	for (i in seq_len(nrow(vertices))) {
		incoming = vertices[i, ] - vertices[previous_index[i], ]
		outgoing = vertices[next_index[i], ] - vertices[i, ]
		incoming_normal = normalize_2d(c(incoming[2], -incoming[1]))
		outgoing_normal = normalize_2d(c(outgoing[2], -outgoing[1]))
		normals[i, ] = normalize_2d(incoming_normal + outgoing_normal)
	}

	return(list(vertices = vertices, normals = normals, v = v))
}

#' @keywords internal
normalize_cross_section_resolution = function(cross_section_resolution) {
	if (length(cross_section_resolution) != 1L || is.na(cross_section_resolution)) {
		stop("cross_section_resolution must be an integer greater than or equal to 8")
	}

	cross_section_resolution = as.integer(cross_section_resolution)
	if (
		!is.finite(cross_section_resolution) ||
			cross_section_resolution < 8L
	) {
		stop("cross_section_resolution must be an integer greater than or equal to 8")
	}

	return(cross_section_resolution)
}

#' @keywords internal
hermite_position = function(p1, p2, m1, m2, s) {
	h00 = 2 * s^3 - 3 * s^2 + 1
	h10 = s^3 - 2 * s^2 + s
	h01 = -2 * s^3 + 3 * s^2
	h11 = s^3 - s^2
	return(h00 * p1 + h10 * m1 + h01 * p2 + h11 * m2)
}

#' @keywords internal
hermite_derivative = function(p1, p2, m1, m2, s) {
	h00 = 6 * s^2 - 6 * s
	h10 = 3 * s^2 - 4 * s + 1
	h01 = -6 * s^2 + 6 * s
	h11 = 3 * s^2 - 2 * s
	return(h00 * p1 + h10 * m1 + h01 * p2 + h11 * m2)
}

#' @keywords internal
parallel_transport_normal = function(
	previous_normal,
	previous_tangent,
	current_tangent
) {
	axis = cross(previous_tangent, current_tangent)
	axis_length = vector_length(axis)
	tangent_dot = clamp_value(sum(previous_tangent * current_tangent), -1, 1)

	if (axis_length <= 1e-8 || tangent_dot >= 1 - 1e-8) {
		projected = project_to_plane(previous_normal, current_tangent)
		if (vector_length(projected) <= 1e-8) {
			return(stable_perpendicular(current_tangent))
		}
		return(normalize_vector(projected, "Failed to transport the ribbon frame"))
	}

	axis = axis / axis_length
	angle = acos(tangent_dot)
	rotated = rodrigues_rotate(previous_normal, axis, angle)
	projected = project_to_plane(rotated, current_tangent)
	if (vector_length(projected) <= 1e-8) {
		return(stable_perpendicular(current_tangent))
	}
	return(normalize_vector(projected, "Failed to transport the ribbon frame"))
}

#' @keywords internal
project_to_plane = function(vector, normal) {
	return(vector - sum(vector * normal) * normal)
}

#' @keywords internal
rodrigues_rotate = function(vector, axis, angle) {
	cosine = cos(angle)
	sine = sin(angle)
	return(
		vector *
			cosine +
			cross(axis, vector) * sine +
			axis * sum(axis * vector) * (1 - cosine)
	)
}

#' @keywords internal
stable_perpendicular = function(vector) {
	axis = c(1, 0, 0)
	if (abs(vector[1]) < abs(vector[2])) {
		axis = c(0, 1, 0)
	}
	if (abs(vector[3]) < min(abs(vector[1]), abs(vector[2]))) {
		axis = c(0, 0, 1)
	}
	perpendicular = cross(vector, axis)
	if (vector_length(perpendicular) <= 1e-8) {
		axis = c(0, 1, 0)
		perpendicular = cross(vector, axis)
	}
	return(normalize_vector(
		perpendicular,
		"Failed to compute a stable ribbon normal"
	))
}

#' @keywords internal
normalize_vector = function(vector, error_message) {
	magnitude = vector_length(vector)
	if (!is.finite(magnitude) || magnitude <= 1e-8) {
		stop(error_message)
	}
	return(vector / magnitude)
}

#' @keywords internal
normalize_2d = function(vector) {
	magnitude = sqrt(sum(vector * vector))
	if (!is.finite(magnitude) || magnitude <= 1e-8) {
		stop("Failed to compute the ribbon cross-section normals")
	}
	return(vector / magnitude)
}

#' @keywords internal
vector_length = function(vector) {
	return(sqrt(sum(vector * vector)))
}

#' @keywords internal
polygon_area = function(vertices) {
	next_vertices = vertices[c(2:nrow(vertices), 1L), , drop = FALSE]
	return(
		sum(
			vertices[, 1] * next_vertices[, 2] - next_vertices[, 1] * vertices[, 2]
		) /
			2
	)
}

#' @keywords internal
clamp_value = function(x, lower, upper) {
	return(max(lower, min(upper, x)))
}

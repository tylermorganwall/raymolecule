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
catmull_rom_chain_point = function(chain, parameter) {
  segment = catmull_rom_parameter_segment(chain, parameter)
  return(hermite_position(
    p1 = segment$segment$p1,
    p2 = segment$segment$p2,
    m1 = segment$segment$m1,
    m2 = segment$segment$m2,
    s = segment$fraction
  ))
}

#' @keywords internal
catmull_rom_chain_derivative = function(chain, parameter) {
  segment = catmull_rom_parameter_segment(chain, parameter)
  return(hermite_derivative(
    p1 = segment$segment$p1,
    p2 = segment$segment$p2,
    m1 = segment$segment$m1,
    m2 = segment$segment$m2,
    s = segment$fraction
  ))
}

#' @keywords internal
catmull_rom_parameter_segment = function(chain, parameter) {
  if (
    length(parameter) != 1L ||
      is.na(parameter) ||
      !is.finite(parameter)
  ) {
    stop("parameter must be a finite number")
  }
  segment_count = length(chain$segments)
  if (segment_count < 1L) {
    stop("catmull_rom_parameter_segment() requires at least one segment")
  }

  parameter = pmax(0, pmin(segment_count, parameter))
  if (parameter >= segment_count) {
    return(list(
      segment = chain$segments[[segment_count]],
      fraction = 1
    ))
  }

  segment_index = floor(parameter) + 1L
  list(
    segment = chain$segments[[segment_index]],
    fraction = parameter - floor(parameter)
  )
}

#' @keywords internal
build_peptide_plane_guides = function(residues) {
  if (nrow(residues) == 0L) {
    stop("build_peptide_plane_guides() requires at least one residue")
  }

  segment_count = nrow(residues) - 1L
  if (segment_count < 1L) {
    return(matrix(NA_real_, nrow = 0L, ncol = 3L))
  }

  segment_guides = matrix(NA_real_, nrow = segment_count, ncol = 3L)
  previous_guide = NULL

  for (segment_index in seq_len(segment_count)) {
    if (
      !residues$has_ca[segment_index] ||
        !residues$has_ca[segment_index + 1L] ||
        !residues$has_o[segment_index]
    ) {
      next
    }

    forward = c(
      residues$ca_x[segment_index + 1L] - residues$ca_x[segment_index],
      residues$ca_y[segment_index + 1L] - residues$ca_y[segment_index],
      residues$ca_z[segment_index + 1L] - residues$ca_z[segment_index]
    )
    guide = c(
      residues$o_x[segment_index] - residues$ca_x[segment_index],
      residues$o_y[segment_index] - residues$ca_y[segment_index],
      residues$o_z[segment_index] - residues$ca_z[segment_index]
    )
    projected = project_to_plane(guide, forward)
    if (vector_length(projected) <= 1e-8) {
      next
    }

    projected = normalize_vector(
      projected,
      error_message = "Failed to construct peptide-plane ribbon guides"
    )
    if (!is.null(previous_guide) && sum(projected * previous_guide) < 0) {
      projected = -projected
    }

    segment_guides[segment_index, ] = projected
    previous_guide = projected
  }

  return(segment_guides)
}

#' @keywords internal
smooth_sheet_control_guides = function(
  control_guides,
  residues,
  kernel = c(1, 4, 6, 4, 1)
) {
  if (
    !is.matrix(control_guides) ||
      ncol(control_guides) != 3L ||
      nrow(control_guides) != nrow(residues)
  ) {
    stop(
      "smooth_sheet_control_guides() requires guides matching the residue table"
    )
  }

  smoothed = control_guides
  sheet_runs = find_sheet_runs(residues)
  if (length(sheet_runs) == 0L) {
    return(smoothed)
  }

  kernel = as.numeric(kernel)
  half_window = as.integer((length(kernel) - 1L) / 2L)

  for (run in sheet_runs) {
    run_rows = run$start:run$end
    valid_rows = run_rows[
      rowSums(is.finite(control_guides[run_rows, , drop = FALSE])) == 3L
    ]
    if (length(valid_rows) < 3L) {
      next
    }

    run_guides = control_guides[valid_rows, , drop = FALSE]
    for (guide_index in 2:nrow(run_guides)) {
      if (sum(run_guides[guide_index, ] * run_guides[guide_index - 1L, ]) < 0) {
        run_guides[guide_index, ] = -run_guides[guide_index, ]
      }
    }

    for (guide_index in seq_len(nrow(run_guides))) {
      window_indices = pmax(
        1L,
        pmin(
          nrow(run_guides),
          (guide_index - half_window):(guide_index + half_window)
        )
      )
      window_guides = run_guides[window_indices, , drop = FALSE]
      averaged = colSums(
        window_guides *
          matrix(kernel, nrow = nrow(window_guides), ncol = 3L)
      ) /
        sum(kernel)
      if (vector_length(averaged) <= 1e-8) {
        next
      }
      smoothed[valid_rows[guide_index], ] = normalize_vector(
        averaged,
        error_message = "Failed to smooth sheet ribbon guides"
      )
    }
  }

  return(smoothed)
}

#' @keywords internal
build_control_guides = function(residues) {
  if (nrow(residues) == 0L) {
    stop("build_control_guides() requires at least one residue")
  }

  control_guides = matrix(NA_real_, nrow = nrow(residues), ncol = 3L)
  segment_guides = build_peptide_plane_guides(residues)

  if (nrow(segment_guides) > 0L) {
    for (residue_index in seq_len(nrow(residues))) {
      candidate_rows = c(residue_index - 1L, residue_index)
      candidate_rows = candidate_rows[
        candidate_rows >= 1L & candidate_rows <= nrow(segment_guides)
      ]
      candidate_guides = segment_guides[candidate_rows, , drop = FALSE]
      valid_candidates = rowSums(is.finite(candidate_guides)) == 3L
      if (!any(valid_candidates)) {
        next
      }
      averaged = colMeans(candidate_guides[valid_candidates, , drop = FALSE])
      if (vector_length(averaged) <= 1e-8) {
        next
      }
      control_guides[residue_index, ] = normalize_vector(
        averaged,
        error_message = "Failed to construct residue ribbon guides"
      )
    }
  }

  fallback_rows = which(
    rowSums(is.finite(control_guides)) != 3L &
      residues$has_ca &
      residues$has_o
  )
  if (length(fallback_rows) > 0L) {
    raw_guides = cbind(
      residues$o_x[fallback_rows] - residues$ca_x[fallback_rows],
      residues$o_y[fallback_rows] - residues$ca_y[fallback_rows],
      residues$o_z[fallback_rows] - residues$ca_z[fallback_rows]
    )
    for (row_index in seq_len(nrow(raw_guides))) {
      control_guides[fallback_rows[row_index], ] = normalize_vector(
        raw_guides[row_index, ],
        error_message = "Failed to construct fallback residue ribbon guides"
      )
    }
  }

  return(smooth_sheet_control_guides(control_guides, residues))
}

#' @keywords internal
build_residue_guides = function(residues, residue_parameter) {
  if (nrow(residues) == 0L) {
    stop("build_residue_guides() requires at least one residue")
  }

  guides = matrix(NA_real_, nrow = length(residue_parameter), ncol = 3L)
  control_guides = build_control_guides(residues)
  valid_guides = rowSums(is.finite(control_guides)) == 3L

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
  residue_parameter = seq_len(nrow(centerline)) - 1,
  residues = NULL
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
    residue_parameter = residue_parameter,
    residues = residues
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
  residue_parameter,
  residues = NULL
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

  knot_indices = knot_indices[valid_knots]
  knot_angles = knot_angles[valid_knots]
  knot_positions = arc_length[knot_indices]

  if (!is.null(residues)) {
    knot_residue_index = pmin(
      nrow(residues),
      pmax(1L, as.integer(round(residue_parameter[knot_indices])) + 1L)
    )
    knot_angles = smooth_sheet_twist_knots(
      knot_positions = knot_positions,
      knot_angles = knot_angles,
      knot_residue_index = knot_residue_index,
      residues = residues
    )
  }

  return(interpolate_twist_angles(
    positions = arc_length,
    knot_positions = knot_positions,
    knot_angles = knot_angles
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

  if (
    abs(flipped_angle - previous_angle) < abs(directed_angle - previous_angle)
  ) {
    return(flipped_angle)
  }

  return(directed_angle)
}

#' @keywords internal
smooth_sheet_twist_knots = function(
  knot_positions,
  knot_angles,
  knot_residue_index,
  residues,
  lambda = 40
) {
  if (
    length(knot_positions) != length(knot_angles) ||
      length(knot_angles) != length(knot_residue_index)
  ) {
    stop("smooth_sheet_twist_knots() requires matching knot vectors")
  }
  if (!is.data.frame(residues)) {
    stop("smooth_sheet_twist_knots() requires a residue data frame")
  }
  if (length(knot_angles) < 3L || nrow(residues) == 0L) {
    return(knot_angles)
  }

  ss_class = residues$ss_class[knot_residue_index]
  sheet_keys = rep("", length(knot_angles))
  sheet_mask = ss_class == "sheet"

  if (!any(sheet_mask)) {
    return(knot_angles)
  }

  sheet_strand = if ("sheet_strand" %in% names(residues)) {
    residues$sheet_strand[knot_residue_index]
  } else {
    rep(NA_integer_, length(knot_angles))
  }
  sheet_keys[sheet_mask] = paste(
    residues$sheet_id[knot_residue_index][sheet_mask],
    sheet_strand[sheet_mask],
    sep = ":"
  )

  smoothed = knot_angles
  index = 1L
  while (index <= length(sheet_keys)) {
    if (!nzchar(sheet_keys[index]) || !is.finite(smoothed[index])) {
      index = index + 1L
      next
    }

    run_end = index
    while (
      run_end < length(sheet_keys) &&
        identical(sheet_keys[run_end + 1L], sheet_keys[index]) &&
        is.finite(smoothed[run_end + 1L])
    ) {
      run_end = run_end + 1L
    }

    run_length = run_end - index + 1L
    if (run_length >= 4L) {
      smoothed[index:run_end] = smooth_twist_angles(
        positions = knot_positions[index:run_end],
        angles = smoothed[index:run_end],
        lambda = lambda
      )
    } else if (run_length == 3L) {
      smoothed[index + 1L] = mean(smoothed[index:run_end])
    }

    index = run_end + 1L
  }

  return(smoothed)
}

#' @keywords internal
interpolate_twist_angles = function(positions, knot_positions, knot_angles) {
  if (
    length(knot_positions) != length(knot_angles) ||
      length(knot_positions) == 0L
  ) {
    stop(
      "interpolate_twist_angles() requires matching non-empty knot positions and angles"
    )
  }
  if (
    any(!is.finite(positions)) ||
      any(!is.finite(knot_positions)) ||
      any(!is.finite(knot_angles))
  ) {
    stop("interpolate_twist_angles() requires finite positions and angles")
  }
  if (any(diff(knot_positions) <= 0)) {
    stop(
      "interpolate_twist_angles() requires strictly increasing knot positions"
    )
  }

  if (length(knot_positions) == 1L) {
    return(rep(knot_angles, length(positions)))
  }

  if (length(knot_positions) == 2L) {
    return(
      stats::approx(
        x = knot_positions,
        y = knot_angles,
        xout = positions,
        method = "linear",
        rule = 2
      )$y
    )
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
      if ("sheet_strand" %in% names(residues)) residues$sheet_strand else
        NA_integer_,
      sep = ":"
    )
  )
  runs = list()
  run_start = NA_integer_
  run_key = NA_character_

  for (i in seq_len(nrow(residues))) {
    if (!sheet_mask[i]) {
      if (!is.na(run_start)) {
        runs[[length(runs) + 1L]] = list(
          start = run_start,
          end = i - 1L,
          key = run_key
        )
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
      runs[[length(runs) + 1L]] = list(
        start = run_start,
        end = i - 1L,
        key = run_key
      )
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
build_peptide_planes = function(residues) {
  if (!is.data.frame(residues)) {
    stop("build_peptide_planes() requires a residue data frame")
  }
  if (nrow(residues) < 2L) {
    stop("build_peptide_planes() requires at least two residues")
  }

  plane_count = nrow(residues) - 1L
  position = matrix(NA_real_, nrow = plane_count, ncol = 3L)
  side = matrix(NA_real_, nrow = plane_count, ncol = 3L)
  normal = matrix(NA_real_, nrow = plane_count, ncol = 3L)
  forward = matrix(NA_real_, nrow = plane_count, ncol = 3L)
  residue1_type = character(plane_count)
  residue2_type = character(plane_count)
  residue3_type = character(plane_count)

  previous_side = NULL
  for (plane_index in seq_len(plane_count)) {
    if (
      !residues$has_ca[plane_index] ||
        !residues$has_ca[plane_index + 1L] ||
        !residues$has_o[plane_index]
    ) {
      stop("Peptide-plane ribbon generation requires contiguous CA and O atoms")
    }

    ca1 = c(
      residues$ca_x[plane_index],
      residues$ca_y[plane_index],
      residues$ca_z[plane_index]
    )
    ca2 = c(
      residues$ca_x[plane_index + 1L],
      residues$ca_y[plane_index + 1L],
      residues$ca_z[plane_index + 1L]
    )
    o1 = c(
      residues$o_x[plane_index],
      residues$o_y[plane_index],
      residues$o_z[plane_index]
    )

    plane_forward = normalize_vector(
      ca2 - ca1,
      error_message = "Failed to construct peptide-plane forward vectors"
    )
    guide = project_to_plane(o1 - ca1, plane_forward)
    if (vector_length(guide) <= 1e-8) {
      stop(
        "Peptide-plane ribbon generation requires non-collinear CA and O atoms"
      )
    }
    plane_side = normalize_vector(
      guide,
      error_message = "Failed to construct peptide-plane side vectors"
    )
    plane_normal = normalize_vector(
      cross(plane_forward, plane_side),
      error_message = "Failed to construct peptide-plane normals"
    )

    if (!is.null(previous_side) && sum(plane_side * previous_side) < 0) {
      plane_side = -plane_side
      plane_normal = -plane_normal
    }

    position[plane_index, ] = (ca1 + ca2) / 2
    side[plane_index, ] = plane_side
    normal[plane_index, ] = plane_normal
    forward[plane_index, ] = plane_forward
    residue1_type[plane_index] = residues$ss_class[plane_index]
    residue2_type[plane_index] = residues$ss_class[min(
      plane_index + 1L,
      nrow(residues)
    )]
    residue3_type[plane_index] = residues$ss_class[min(
      plane_index + 2L,
      nrow(residues)
    )]
    previous_side = plane_side
  }

  return(list(
    position = position,
    side = side,
    normal = normal,
    forward = forward,
    residue1_type = residue1_type,
    residue2_type = residue2_type,
    residue3_type = residue3_type
  ))
}

#' @keywords internal
ss_class_to_ribbon_type = function(ss_class) {
  if (identical(ss_class, "sheet")) {
    return("strand")
  }
  if (identical(ss_class, "helix")) {
    return("helix")
  }
  return("coil")
}

#' @keywords internal
transition_ribbon_types = function(planes, plane_index) {
  type1 = ss_class_to_ribbon_type(planes$residue2_type[plane_index])
  type2 = type1
  type0 = ss_class_to_ribbon_type(planes$residue1_type[plane_index])
  type3 = ss_class_to_ribbon_type(planes$residue3_type[plane_index])
  type_rank = c(coil = 1L, helix = 2L, strand = 3L)

  if (type_rank[[type1]] > type_rank[[type0]] && identical(type1, type3)) {
    type1 = type0
  }
  if (
    type_rank[[ss_class_to_ribbon_type(planes$residue2_type[plane_index])]] >
      type_rank[[type3]] &&
      identical(
        type0,
        ss_class_to_ribbon_type(planes$residue2_type[plane_index])
      )
  ) {
    type2 = type3
  }

  return(list(type0 = type0, type1 = type1, type2 = type2))
}

#' @keywords internal
ellipse_profile_2d = function(n, width, height) {
  theta = seq(0, 2 * pi, length.out = n + 1L)
  theta = theta[-length(theta)] + pi / 4
  return(canonicalize_profile_phase(cbind(
    cos(theta) * width / 2,
    sin(theta) * height / 2
  )))
}

#' @keywords internal
rectangle_profile_2d = function(n, width, height) {
  half_width = width / 2
  half_height = height / 2
  counts = rep(floor(n / 4), 4L)
  remainder = n - sum(counts)
  if (remainder > 0L) {
    counts[seq_len(remainder)] = counts[seq_len(remainder)] + 1L
  }

  corners = rbind(
    c(half_width, half_height),
    c(-half_width, half_height),
    c(-half_width, -half_height),
    c(half_width, -half_height)
  )
  vertices = matrix(NA_real_, nrow = n, ncol = 2L)
  vertex_index = 1L
  for (segment_index in seq_len(4L)) {
    count = counts[segment_index]
    start = corners[segment_index, ]
    end = corners[if (segment_index == 4L) 1L else segment_index + 1L, ]
    t_values = seq(0, 1, length.out = count + 1L)[seq_len(count)]
    for (t in t_values) {
      vertices[vertex_index, ] = (1 - t) * start + t * end
      vertex_index = vertex_index + 1L
    }
  }
  return(canonicalize_profile_phase(vertices))
}

#' @keywords internal
rounded_rectangle_profile_2d = function(n, width, height) {
  return(canonicalize_profile_phase(
    ribbon_cross_section_profile(
      ribbon_width = width,
      ribbon_thickness = height,
      cross_section_resolution = n,
      shape_exponent = 4
    )$vertices
  ))
}

#' @keywords internal
canonicalize_profile_phase = function(profile_vertices) {
  if (!is.matrix(profile_vertices) || ncol(profile_vertices) != 2L) {
    stop("canonicalize_profile_phase() requires an n x 2 profile matrix")
  }

  start_scores = profile_vertices[, 1] + profile_vertices[, 2]
  start_index = which.max(start_scores)
  return(profile_vertices[
    c(start_index:nrow(profile_vertices), seq_len(start_index - 1L)),
    ,
    drop = FALSE
  ])
}

#' @keywords internal
ribbon_segment_profiles = function(
  planes,
  plane_index,
  cross_section_resolution,
  ribbon_width,
  ribbon_thickness
) {
  types = transition_ribbon_types(planes, plane_index)
  type0 = types$type0
  type1 = types$type1
  type2 = types$type2

  arrow_body_width = ribbon_width
  arrow_head_width = ribbon_width * 1.5
  arrow_height = max(ribbon_thickness * 1.75, ribbon_width * 0.18)
  arrow_tip_width = max(ribbon_width * 0.02, 1e-3)
  tube_size = max(ribbon_thickness * 2.25, ribbon_width * 0.32)

  make_profile = function(type, width_override = NULL) {
    if (identical(type, "strand")) {
      return(rectangle_profile_2d(
        cross_section_resolution,
        if (is.null(width_override)) arrow_body_width else width_override,
        arrow_height
      ))
    }
    if (identical(type, "helix")) {
      return(rounded_rectangle_profile_2d(
        cross_section_resolution,
        ribbon_width,
        ribbon_thickness
      ))
    }
    return(ellipse_profile_2d(cross_section_resolution, tube_size, tube_size))
  }
  body_profile = make_profile("strand")
  head_profile = make_profile("strand", width_override = arrow_head_width)
  tip_profile = make_profile("strand", width_override = arrow_tip_width)

  if (identical(type0, "strand") && !identical(type1, "strand")) {
    profile1 = tip_profile
  } else if (identical(type1, "strand")) {
    if (identical(type2, "strand")) {
      profile1 = body_profile
    } else {
      profile1 = body_profile
    }
  } else if (identical(type1, "helix")) {
    profile1 = make_profile("helix")
  } else {
    profile1 = make_profile("coil")
  }

  if (identical(type2, "strand")) {
    profile2 = body_profile
  } else if (identical(type2, "helix")) {
    profile2 = make_profile("helix")
  } else {
    profile2 = make_profile("coil")
  }

  if (identical(type1, "strand") && !identical(type2, "strand")) {
    profile2 = tip_profile
  }

  return(list(
    start = profile1,
    end = profile2,
    body = body_profile,
    head = head_profile,
    tip = tip_profile,
    type0 = type0,
    type1 = type1,
    type2 = type2
  ))
}

#' @keywords internal
segment_profile_fraction = function(type0, type1, type2, fraction) {
  if (identical(type1, "strand") && !identical(type2, "strand")) {
    return(fraction)
  }
  if (identical(type0, "strand") && !identical(type1, "strand")) {
    return(sqrt(1 - (1 - fraction)^2))
  }
  return(smoothstep01(fraction))
}

#' @keywords internal
interpolate_segment_profile = function(profiles, fraction) {
  if (
    identical(profiles$type1, "strand") && !identical(profiles$type2, "strand")
  ) {
    base_step = strand_arrow_base_step()

    if (fraction < base_step) {
      return(profiles$body)
    }
    u = (fraction - base_step) / (1 - base_step)
    u = pmax(0, pmin(1, u))
    return((1 - u) * profiles$head + u * profiles$tip)
  }

  profile_fraction = segment_profile_fraction(
    type0 = profiles$type0,
    type1 = profiles$type1,
    type2 = profiles$type2,
    fraction = fraction
  )
  return(
    (1 - profile_fraction) * profiles$start + profile_fraction * profiles$end
  )
}

#' @keywords internal
segment_helix_centerline_fraction = function(profiles, fraction) {
  start_helix = identical(profiles$type1, "helix")
  end_helix = identical(profiles$type2, "helix")

  if (start_helix && end_helix) {
    return(1)
  }
  if (!start_helix && !end_helix) {
    return(0)
  }

  transition = smoothstep01(fraction)
  if (start_helix) {
    return(1 - transition)
  }
  return(transition)
}

#' @keywords internal
segment_helix_centerline_fraction_derivative = function(profiles, fraction) {
  start_helix = identical(profiles$type1, "helix")
  end_helix = identical(profiles$type2, "helix")

  if (start_helix == end_helix) {
    return(0)
  }

  clamped = pmax(0, pmin(1, fraction))
  transition_derivative = 6 * clamped * (1 - clamped)
  if (start_helix) {
    return(-transition_derivative)
  }
  return(transition_derivative)
}

#' @keywords internal
build_helix_centerline_path = function(residues) {
  helix_rows = which(residues$ss_class == "helix" & residues$has_ca)
  if (length(helix_rows) < 5L) {
    return(list(runs = list()))
  }

  run_breaks = c(TRUE, diff(helix_rows) != 1L)
  run_ids = cumsum(run_breaks)
  run_rows = split(helix_rows, run_ids)
  runs = list()

  for (rows in run_rows) {
    if (length(rows) < 5L) {
      next
    }
    points = cbind(
      residues$ca_x[rows],
      residues$ca_y[rows],
      residues$ca_z[rows]
    )
    run = build_single_helix_centerline_run(
      points = points,
      parameter = rows - 1L
    )
    if (!is.null(run)) {
      runs[[length(runs) + 1L]] = run
    }
  }

  return(list(runs = runs))
}

#' @keywords internal
build_single_helix_centerline_run = function(points, parameter) {
  if (nrow(points) < 5L) {
    return(NULL)
  }

  origin = colMeans(points)
  centered = sweep(points, 2, origin, FUN = "-")
  axis = estimate_helix_axis(points)
  if (is.null(axis)) {
    return(NULL)
  }
  if (sum(axis * (points[nrow(points), ] - points[1, ])) < 0) {
    axis = -axis
  }

  height = drop(centered %*% axis)
  projected = centered -
    height *
      matrix(
        rep(axis, each = nrow(centered)),
        ncol = 3L
      )
  projection_lengths = sqrt(rowSums(projected^2))
  u_axis = projected[which.max(projection_lengths), ]
  if (vector_length(u_axis) <= 1e-8) {
    u_axis = stable_perpendicular(axis)
  }
  u_axis = normalize_vector(
    u_axis,
    error_message = "Failed to construct helix centerline basis"
  )
  v_axis = normalize_vector(
    cross(axis, u_axis),
    error_message = "Failed to construct helix centerline basis"
  )

  x = drop(centered %*% u_axis)
  y = drop(centered %*% v_axis)
  circle = fit_helix_circle_2d(x, y)
  if (is.null(circle)) {
    return(NULL)
  }

  radial_x = x - circle$center[[1]]
  radial_y = y - circle$center[[2]]
  radius = sqrt(radial_x^2 + radial_y^2)
  if (any(!is.finite(radius)) || mean(radius) <= 1e-8) {
    return(NULL)
  }

  angle = unwrap_ribbon_angles(atan2(radial_y, radial_x))
  if (
    !all(is.finite(angle)) || abs(angle[length(angle)] - angle[[1]]) <= 1e-8
  ) {
    return(NULL)
  }
  angle_spline = make_helix_scalar_spline(parameter, angle)
  height_spline = make_helix_scalar_spline(parameter, height)

  return(list(
    parameter = parameter,
    origin = origin,
    axis = axis,
    u_axis = u_axis,
    v_axis = v_axis,
    circle_center = circle$center,
    radius = mean(radius),
    angle = angle,
    height = height,
    angle_spline = angle_spline,
    height_spline = height_spline
  ))
}

#' @keywords internal
make_helix_scalar_spline = function(parameter, values) {
  method = if (
    all(diff(values) >= 0) ||
      all(diff(values) <= 0)
  ) {
    "monoH.FC"
  } else {
    "natural"
  }
  spline = stats::splinefun(parameter, values, method = method)

  list(
    spline = spline,
    first_parameter = parameter[[1]],
    last_parameter = parameter[[length(parameter)]],
    first_value = values[[1]],
    last_value = values[[length(values)]],
    first_slope = (values[[2]] - values[[1]]) /
      (parameter[[2]] - parameter[[1]]),
    last_slope = (values[[length(values)]] - values[[length(values) - 1L]]) /
      (parameter[[length(parameter)]] - parameter[[length(parameter) - 1L]])
  )
}

#' @keywords internal
evaluate_helix_scalar_spline = function(spline, parameter, derivative = 0L) {
  if (parameter < spline$first_parameter) {
    if (derivative == 0L) {
      return(
        spline$first_value +
          (parameter - spline$first_parameter) * spline$first_slope
      )
    }
    return(spline$first_slope)
  }
  if (parameter > spline$last_parameter) {
    if (derivative == 0L) {
      return(
        spline$last_value +
          (parameter - spline$last_parameter) * spline$last_slope
      )
    }
    return(spline$last_slope)
  }

  spline$spline(parameter, deriv = derivative)
}

#' @keywords internal
estimate_helix_axis = function(points) {
  if (nrow(points) < 2L) {
    return(NULL)
  }

  candidate_lags = unique(c(
    seq_len(min(8L, nrow(points) - 1L)),
    nrow(points) - 1L
  ))
  best_axis = NULL
  best_score = Inf

  for (lag in candidate_lags) {
    lag_vectors = points[(lag + 1L):nrow(points), , drop = FALSE] -
      points[seq_len(nrow(points) - lag), , drop = FALSE]
    axis = colSums(lag_vectors)
    if (vector_length(axis) <= 1e-8) {
      next
    }
    axis = normalize_vector(
      axis,
      error_message = "Failed to construct helix centerline axis"
    )
    score = score_helix_axis(points, axis)
    if (is.finite(score) && score < best_score) {
      best_axis = axis
      best_score = score
    }
  }

  if (!is.null(best_axis)) {
    return(best_axis)
  }

  centered = sweep(points, 2, colMeans(points), FUN = "-")
  fit = tryCatch(
    svd(centered),
    error = function(e) NULL
  )
  if (is.null(fit) || ncol(fit$v) < 1L) {
    return(NULL)
  }

  return(normalize_vector(
    fit$v[, 1],
    error_message = "Failed to construct helix centerline axis"
  ))
}

#' @keywords internal
score_helix_axis = function(points, axis) {
  origin = colMeans(points)
  centered = sweep(points, 2, origin, FUN = "-")
  height = drop(centered %*% axis)
  projected = centered -
    height *
      matrix(
        rep(axis, each = nrow(centered)),
        ncol = 3L
      )
  projection_lengths = sqrt(rowSums(projected^2))
  if (!any(is.finite(projection_lengths) & projection_lengths > 1e-8)) {
    return(Inf)
  }

  u_axis = projected[which.max(projection_lengths), ]
  u_axis = normalize_vector(
    u_axis,
    error_message = "Failed to construct helix centerline basis"
  )
  v_axis = normalize_vector(
    cross(axis, u_axis),
    error_message = "Failed to construct helix centerline basis"
  )
  x = drop(centered %*% u_axis)
  y = drop(centered %*% v_axis)
  circle = fit_helix_circle_2d(x, y)
  if (is.null(circle)) {
    return(Inf)
  }

  radius = sqrt((x - circle$center[[1]])^2 + (y - circle$center[[2]])^2)
  if (any(!is.finite(radius)) || mean(radius) <= 1e-8) {
    return(Inf)
  }

  sd(radius) / mean(radius)
}

#' @keywords internal
fit_helix_circle_2d = function(x, y) {
  design = cbind(x, y, 1)
  fit = tryCatch(
    stats::lm.fit(design, -(x^2 + y^2)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$rank < 3L || any(!is.finite(fit$coefficients))) {
    center = c(mean(x), mean(y))
    radius = mean(sqrt((x - center[[1]])^2 + (y - center[[2]])^2))
    if (!is.finite(radius) || radius <= 1e-8) {
      return(NULL)
    }
    return(list(center = center, radius = radius))
  }

  coefficients = fit$coefficients
  center = -coefficients[1:2] / 2
  radius_squared = sum(center^2) - coefficients[[3]]
  if (!is.finite(radius_squared) || radius_squared <= 1e-12) {
    return(NULL)
  }

  return(list(center = center, radius = sqrt(radius_squared)))
}

#' @keywords internal
unwrap_ribbon_angles = function(angles) {
  if (length(angles) < 2L) {
    return(angles)
  }

  unwrapped = angles
  for (i in 2:length(unwrapped)) {
    while (unwrapped[[i]] - unwrapped[[i - 1L]] > pi) {
      unwrapped[[i]] = unwrapped[[i]] - 2 * pi
    }
    while (unwrapped[[i]] - unwrapped[[i - 1L]] < -pi) {
      unwrapped[[i]] = unwrapped[[i]] + 2 * pi
    }
  }
  return(unwrapped)
}

#' @keywords internal
sample_helix_centerline_path = function(path, parameter) {
  if (length(path$runs) == 0L) {
    return(NULL)
  }

  for (run in path$runs) {
    first_parameter = run$parameter[[1]]
    last_parameter = run$parameter[[length(run$parameter)]]
    entry_transition_margin = 1.5
    exit_transition_margin = 0.5
    if (
      parameter < first_parameter - entry_transition_margin ||
        parameter > last_parameter + exit_transition_margin
    ) {
      next
    }

    angle = evaluate_helix_scalar_spline(
      spline = run$angle_spline,
      parameter = parameter
    )
    height = evaluate_helix_scalar_spline(
      spline = run$height_spline,
      parameter = parameter
    )
    angle_derivative = evaluate_helix_scalar_spline(
      spline = run$angle_spline,
      parameter = parameter,
      derivative = 1L
    )
    height_derivative = evaluate_helix_scalar_spline(
      spline = run$height_spline,
      parameter = parameter,
      derivative = 1L
    )
    radial = cos(angle) * run$u_axis + sin(angle) * run$v_axis
    angular_direction = -sin(angle) * run$u_axis + cos(angle) * run$v_axis
    center =
      run$origin +
      height * run$axis +
      run$circle_center[[1]] * run$u_axis +
      run$circle_center[[2]] * run$v_axis +
      run$radius * radial
    derivative = run$radius *
      angle_derivative *
      angular_direction +
      height_derivative * run$axis
    if (vector_length(derivative) <= 1e-8) {
      return(NULL)
    }

    return(list(center = center, derivative = derivative))
  }

  return(NULL)
}

#' @keywords internal
strand_arrow_base_step = function() {
  return(0.001)
}

#' @keywords internal
segment_sample_fractions = function(
  profiles,
  subdivisions,
  include_endpoint = TRUE
) {
  fractions = seq(0, 1, length.out = subdivisions + 1L)
  if (!include_endpoint) {
    fractions = fractions[-length(fractions)]
  }
  if (
    identical(profiles$type1, "strand") && !identical(profiles$type2, "strand")
  ) {
    fractions = sort(unique(c(fractions, strand_arrow_base_step())))
  }
  return(fractions)
}

#' @keywords internal
bspline_basis = function(t) {
  return(c(
    (-t^3 + 3 * t^2 - 3 * t + 1) / 6,
    (3 * t^3 - 6 * t^2 + 4) / 6,
    (-3 * t^3 + 3 * t^2 + 3 * t + 1) / 6,
    t^3 / 6
  ))
}

#' @keywords internal
bspline_basis_derivative = function(t) {
  return(c(
    (-3 * t^2 + 6 * t - 3) / 6,
    (9 * t^2 - 12 * t) / 6,
    (-9 * t^2 + 6 * t + 3) / 6,
    (3 * t^2) / 6
  ))
}

#' @keywords internal
bspline_point = function(control_points, fraction) {
  basis = bspline_basis(fraction)
  return(drop(basis %*% control_points))
}

#' @keywords internal
bspline_derivative_point = function(control_points, fraction) {
  basis = bspline_basis_derivative(fraction)
  return(drop(basis %*% control_points))
}

#' @keywords internal
profile_perimeter_fraction = function(profile_vertices) {
  edge_vectors = profile_vertices[
    c(2:nrow(profile_vertices), 1L),
    ,
    drop = FALSE
  ] -
    profile_vertices
  edge_lengths = sqrt(rowSums(edge_vectors^2))
  perimeter = sum(edge_lengths)
  if (!is.finite(perimeter) || perimeter <= 1e-8) {
    return(seq(0, 1, length.out = nrow(profile_vertices) + 1L)[seq_len(nrow(
      profile_vertices
    ))])
  }
  cumulative = c(0, cumsum(edge_lengths))
  return(cumulative[seq_len(nrow(profile_vertices))] / perimeter)
}

#' @keywords internal
profile_vertex_normals_2d = function(profile_vertices) {
  previous_index = c(
    nrow(profile_vertices),
    seq_len(nrow(profile_vertices) - 1L)
  )
  next_index = c(2:nrow(profile_vertices), 1L)
  normals = matrix(NA_real_, nrow = nrow(profile_vertices), ncol = 2L)

  for (i in seq_len(nrow(profile_vertices))) {
    incoming = profile_vertices[i, ] - profile_vertices[previous_index[i], ]
    outgoing = profile_vertices[next_index[i], ] - profile_vertices[i, ]
    incoming_length = sqrt(sum(incoming * incoming))
    outgoing_length = sqrt(sum(outgoing * outgoing))

    if (incoming_length <= 1e-8 && outgoing_length <= 1e-8) {
      reference = profile_vertices[i, ]
      if (sqrt(sum(reference * reference)) <= 1e-8) {
        normals[i, ] = c(1, 0)
      } else {
        normals[i, ] = reference / sqrt(sum(reference * reference))
      }
      next
    }

    incoming_normal = if (incoming_length <= 1e-8) {
      c(0, 0)
    } else {
      c(incoming[2], -incoming[1]) / incoming_length
    }
    outgoing_normal = if (outgoing_length <= 1e-8) {
      c(0, 0)
    } else {
      c(outgoing[2], -outgoing[1]) / outgoing_length
    }
    combined = incoming_normal + outgoing_normal
    if (sqrt(sum(combined * combined)) <= 1e-8) {
      combined = if (outgoing_length > 1e-8) outgoing_normal else
        incoming_normal
    }
    normals[i, ] = combined / sqrt(sum(combined * combined))
  }

  return(normals)
}

#' @keywords internal
build_peptide_surface_chain_mesh = function(
  residues,
  ribbon_width,
  ribbon_thickness,
  cross_section_resolution,
  subdivisions,
  chain_id,
  chain_index
) {
  planes = build_peptide_planes(residues)
  plane_count = nrow(planes$position)
  if (plane_count < 2L) {
    stop(
      "build_peptide_surface_chain_mesh() requires at least two peptide planes"
    )
  }

  ca_chain = catmull_rom_chain(cbind(
    residues$ca_x,
    residues$ca_y,
    residues$ca_z
  ))
  helix_path = build_helix_centerline_path(residues)
  segment_count = plane_count - 1L
  segment_profiles = vector(mode = "list", length = segment_count)
  segment_fraction_list = vector(mode = "list", length = segment_count)
  for (segment_index in seq_len(segment_count)) {
    segment_profiles[[segment_index]] = ribbon_segment_profiles(
      planes = planes,
      plane_index = segment_index,
      cross_section_resolution = cross_section_resolution,
      ribbon_width = ribbon_width,
      ribbon_thickness = ribbon_thickness
    )
    segment_fraction_list[[segment_index]] = segment_sample_fractions(
      profiles = segment_profiles[[segment_index]],
      subdivisions = subdivisions,
      include_endpoint = segment_index == segment_count
    )
  }
  ring_count = sum(vapply(segment_fraction_list, length, integer(1)))
  ring_size = cross_section_resolution
  side_vertex_count = ring_count * ring_size

  vertices = matrix(NA_real_, nrow = side_vertex_count + 2L, ncol = 3L)
  normals = matrix(NA_real_, nrow = side_vertex_count + 2L, ncol = 3L)
  texcoords = matrix(
    NA_real_,
    nrow = ring_count * (ring_size + 1L) + 2L,
    ncol = 2L
  )
  centerline = matrix(NA_real_, nrow = ring_count, ncol = 3L)
  tangent = matrix(NA_real_, nrow = ring_count, ncol = 3L)
  normal_axis = matrix(NA_real_, nrow = ring_count, ncol = 3L)
  binormal_axis = matrix(NA_real_, nrow = ring_count, ncol = 3L)
  residue_parameter = numeric(ring_count)
  sampled_width = numeric(ring_count)
  sampled_thickness = numeric(ring_count)
  profile_kind = character(ring_count)
  ring_v = vector(mode = "list", length = ring_count)

  previous_side = NULL
  ring_counter = 1L
  for (segment_index in seq_len(segment_count)) {
    control_rows = pmax(
      1L,
      pmin(
        plane_count,
        c(
          segment_index - 1L,
          segment_index,
          segment_index + 1L,
          segment_index + 2L
        )
      )
    )
    position_controls = planes$position[control_rows, , drop = FALSE]
    side_controls = planes$side[control_rows, , drop = FALSE]
    normal_controls = planes$normal[control_rows, , drop = FALSE]
    profiles = segment_profiles[[segment_index]]
    s_values = segment_fraction_list[[segment_index]]

    for (fraction in s_values) {
      profile_vertices = interpolate_segment_profile(
        profiles = profiles,
        fraction = fraction
      )
      profile_normals = profile_vertex_normals_2d(profile_vertices)
      residue_parameter_value = segment_index - 0.5 + fraction
      peptide_center = bspline_point(position_controls, fraction)
      peptide_derivative = bspline_derivative_point(position_controls, fraction)
      ca_center = catmull_rom_chain_point(
        chain = ca_chain,
        parameter = residue_parameter_value
      )
      ca_derivative = catmull_rom_chain_derivative(
        chain = ca_chain,
        parameter = residue_parameter_value
      )
      helix_centerline_fraction = segment_helix_centerline_fraction(
        profiles = profiles,
        fraction = fraction
      )
      helix_centerline_fraction_derivative =
        segment_helix_centerline_fraction_derivative(
          profiles = profiles,
          fraction = fraction
        )
      if (helix_centerline_fraction > 0) {
        helix_sample = sample_helix_centerline_path(
          path = helix_path,
          parameter = residue_parameter_value
        )
        if (!is.null(helix_sample)) {
          ca_center = helix_sample$center
          ca_derivative = helix_sample$derivative
        }
      }
      current_center =
        (1 - helix_centerline_fraction) *
        peptide_center +
        helix_centerline_fraction * ca_center
      current_tangent = normalize_vector(
        (1 - helix_centerline_fraction) *
          peptide_derivative +
          helix_centerline_fraction * ca_derivative +
          helix_centerline_fraction_derivative *
            (ca_center - peptide_center),
        error_message = "Peptide-plane ribbon sampling produced a zero tangent vector"
      )
      side_seed = project_to_plane(
        bspline_point(side_controls, fraction),
        current_tangent
      )
      if (vector_length(side_seed) <= 1e-8) {
        side_seed = stable_perpendicular(current_tangent)
      }
      current_side = normalize_vector(
        side_seed,
        error_message = "Failed to construct peptide-plane ribbon side vectors"
      )
      if (!is.null(previous_side) && sum(current_side * previous_side) < 0) {
        current_side = -current_side
      }

      normal_seed = project_to_plane(
        bspline_point(normal_controls, fraction),
        current_tangent
      )
      current_binormal = cross(current_tangent, current_side)
      if (vector_length(current_binormal) <= 1e-8) {
        stop("Failed to construct peptide-plane ribbon normals")
      }
      current_binormal = normalize_vector(
        current_binormal,
        error_message = "Failed to construct peptide-plane ribbon normals"
      )
      if (
        vector_length(normal_seed) > 1e-8 &&
          sum(normal_seed * current_binormal) < 0
      ) {
        current_binormal = -current_binormal
      }
      current_side = normalize_vector(
        cross(current_binormal, current_tangent),
        error_message = "Failed to orthogonalize peptide-plane ribbon axes"
      )
      previous_side = current_side

      vertex_rows = ((ring_counter - 1L) * ring_size + 1L):(ring_counter *
        ring_size)
      ring_vertices = matrix(NA_real_, nrow = ring_size, ncol = 3L)
      ring_vertices =
        matrix(
          rep(current_center, each = ring_size),
          ncol = 3L,
          byrow = FALSE
        ) +
        profile_vertices[, 1] *
          matrix(rep(current_side, each = ring_size), ncol = 3L) +
        profile_vertices[, 2] *
          matrix(rep(current_binormal, each = ring_size), ncol = 3L)

      vertices[vertex_rows, ] = ring_vertices
      normals[vertex_rows, ] =
        profile_normals[, 1] *
        matrix(rep(current_side, each = ring_size), ncol = 3L) +
        profile_normals[, 2] *
          matrix(rep(current_binormal, each = ring_size), ncol = 3L)
      centerline[ring_counter, ] = current_center
      tangent[ring_counter, ] = current_tangent
      normal_axis[ring_counter, ] = current_side
      binormal_axis[ring_counter, ] = current_binormal
      residue_parameter[ring_counter] = residue_parameter_value
      sampled_width[ring_counter] = diff(range(profile_vertices[, 1]))
      sampled_thickness[ring_counter] = diff(range(profile_vertices[, 2]))
      profile_kind[ring_counter] = profiles$type1
      ring_v[[ring_counter]] = profile_perimeter_fraction(profile_vertices)
      ring_counter = ring_counter + 1L
    }
  }

  arc_length = numeric(ring_count)
  if (ring_count > 1L) {
    arc_length[-1] = cumsum(sqrt(rowSums(
      (centerline[-1, , drop = FALSE] -
        centerline[-ring_count, , drop = FALSE])^2
    )))
  }
  if (arc_length[ring_count] <= 0) {
    u_values = rep(0, ring_count)
  } else {
    u_values = arc_length / arc_length[ring_count]
    u_values = pmin(u_values, ribbon_texture_u_max())
  }

  for (ring_index in seq_len(ring_count)) {
    tex_rows = ((ring_index - 1L) * (ring_size + 1L) + 1L):(ring_index *
      (ring_size + 1L))
    texcoords[tex_rows[-length(tex_rows)], ] = cbind(
      rep(u_values[ring_index], ring_size),
      ring_v[[ring_index]]
    )
    texcoords[tex_rows[length(tex_rows)], ] = c(u_values[ring_index], 1)
  }

  normals[seq_len(side_vertex_count), ] = compute_ribbon_surface_normals(
    vertices = vertices[seq_len(side_vertex_count), , drop = FALSE],
    centerline = centerline,
    ring_count = ring_count,
    ring_size = ring_size,
    fallback_normals = normals[seq_len(side_vertex_count), , drop = FALSE]
  )

  start_center_index = side_vertex_count + 1L
  end_center_index = side_vertex_count + 2L
  vertices[start_center_index, ] = centerline[1, ]
  vertices[end_center_index, ] = centerline[ring_count, ]
  normals[start_center_index, ] = -tangent[1, ]
  normals[end_center_index, ] = tangent[ring_count, ]

  start_center_tex = ring_count * (ring_size + 1L) + 1L
  end_center_tex = start_center_tex + 1L
  texcoords[start_center_tex, ] = c(0, 0.5)
  texcoords[end_center_tex, ] = c(ribbon_texture_u_max(), 0.5)

  side_data = connect_ribbon_rings(
    vertices = vertices[seq_len(side_vertex_count), , drop = FALSE],
    ring_count = ring_count,
    ring_size = ring_size,
    normals = normals[seq_len(side_vertex_count), , drop = FALSE]
  )
  cap_data = cap_ribbon_ends(ring_count = ring_count, ring_size = ring_size)

  sampled_df = data.frame(
    x = centerline[, 1],
    y = centerline[, 2],
    z = centerline[, 3],
    tangent_x = tangent[, 1],
    tangent_y = tangent[, 2],
    tangent_z = tangent[, 3],
    normal_x = normal_axis[, 1],
    normal_y = normal_axis[, 2],
    normal_z = normal_axis[, 3],
    binormal_x = binormal_axis[, 1],
    binormal_y = binormal_axis[, 2],
    binormal_z = binormal_axis[, 3],
    arc_length = arc_length,
    chain_id = rep(chain_id, ring_count),
    residue_parameter = residue_parameter,
    ribbon_width = sampled_width,
    ribbon_thickness = sampled_thickness,
    profile_kind = profile_kind,
    stringsAsFactors = FALSE
  )

  return(list(
    chain_id = chain_id,
    chain_index = chain_index,
    vertices = vertices,
    indices = rbind(side_data$indices, cap_data$indices),
    normals = normals,
    norm_indices = rbind(side_data$norm_indices, cap_data$norm_indices),
    texcoords = texcoords,
    tex_indices = rbind(side_data$tex_indices, cap_data$tex_indices),
    sampled = list(
      centerline = centerline,
      tangent = tangent,
      normal = normal_axis,
      binormal = binormal_axis,
      arc_length = arc_length,
      chain_id = rep(chain_id, ring_count),
      residue_parameter = residue_parameter,
      sampled = sampled_df
    )
  ))
}

#' @keywords internal
peptide_surface_segment_eligible = function(residues) {
  if (!is.data.frame(residues)) {
    stop("peptide_surface_segment_eligible() requires a residue data frame")
  }
  if (nrow(residues) < 3L) {
    return(FALSE)
  }
  if (!all(residues$has_ca)) {
    return(FALSE)
  }
  if (!all(residues$has_o[seq_len(nrow(residues) - 1L)])) {
    return(FALSE)
  }
  return(TRUE)
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
    use_peptide_surface = peptide_surface_segment_eligible(chain_residues)

    if (use_peptide_surface) {
      chain_meshes[[segment_index]] = build_peptide_surface_chain_mesh(
        residues = chain_residues,
        ribbon_width = ribbon_width,
        ribbon_thickness = ribbon_thickness,
        cross_section_resolution = cross_section_resolution,
        subdivisions = subdivisions,
        chain_id = chain_residues$chain_id[1],
        chain_index = segment_index
      )
      sampled_paths[[segment_index]] = chain_meshes[[segment_index]]$sampled
      next
    }

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
      residue_parameter = sampled_chain$residue_parameter,
      residues = chain_residues
    )
    ribbon_widths = rep(ribbon_width, length(frame$residue_parameter))
    ribbon_thicknesses = rep(ribbon_thickness, length(frame$residue_parameter))
    shape_exponents = rep(4, length(frame$residue_parameter))
    frame$sampled$ribbon_width = ribbon_widths
    frame$sampled$ribbon_thickness = ribbon_thicknesses
    frame$sampled$ribbon_profile_exponent = shape_exponents

    chain_meshes[[segment_index]] = build_single_ribbon_chain_mesh(
      frame = frame,
      ribbon_widths = ribbon_widths,
      ribbon_thicknesses = ribbon_thicknesses,
      profile_exponents = shape_exponents,
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
  if (
    !is.matrix(centerline) ||
      ncol(centerline) != 3L ||
      nrow(centerline) != ring_count
  ) {
    stop(
      "compute_ribbon_surface_normals() requires centerline rows matching ring_count"
    )
  }
  if (nrow(vertices) != ring_count * ring_size) {
    stop(
      "compute_ribbon_surface_normals() requires vertices matching ring_count * ring_size"
    )
  }
  if (!is.null(fallback_normals)) {
    if (
      !is.matrix(fallback_normals) ||
        ncol(fallback_normals) != 3L ||
        nrow(fallback_normals) != nrow(vertices)
    ) {
      stop(
        "compute_ribbon_surface_normals() fallback_normals must match vertices"
      )
    }
  }

  normals = matrix(NA_real_, nrow = nrow(vertices), ncol = 3L)

  for (ring_index in seq_len(ring_count)) {
    prev_ring = max(1L, ring_index - 1L)
    next_ring = min(ring_count, ring_index + 1L)

    for (vertex_index in seq_len(ring_size)) {
      current_row = (ring_index - 1L) * ring_size + vertex_index
      prev_vertex_index = if (vertex_index == 1L) ring_size else
        vertex_index - 1L
      next_vertex_index = if (vertex_index == ring_size) 1L else
        vertex_index + 1L
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

connect_ribbon_rings = function(
  vertices,
  ring_count,
  ring_size,
  normals = NULL
) {
  if (!is.null(normals)) {
    if (
      !is.matrix(normals) ||
        ncol(normals) != 3L ||
        nrow(normals) != nrow(vertices)
    ) {
      stop("connect_ribbon_rings() normals must match the vertex matrix")
    }
  }

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

      if (is.null(normals)) {
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
      } else {
        candidate_indices = list(
          rbind(c(a, b, c), c(a, c, d)),
          rbind(c(a, c, b), c(a, d, c)),
          rbind(c(a, b, d), c(b, c, d)),
          rbind(c(a, d, b), c(b, d, c))
        )
        candidate_tex = list(
          rbind(
            c(tex_a + vertex_index, b_tex, c_tex),
            c(tex_a + vertex_index, c_tex, tex_b + vertex_index)
          ),
          rbind(
            c(tex_a + vertex_index, c_tex, b_tex),
            c(tex_a + vertex_index, tex_b + vertex_index, c_tex)
          ),
          rbind(
            c(tex_a + vertex_index, b_tex, tex_b + vertex_index),
            c(b_tex, c_tex, tex_b + vertex_index)
          ),
          rbind(
            c(tex_a + vertex_index, tex_b + vertex_index, b_tex),
            c(b_tex, tex_b + vertex_index, c_tex)
          )
        )
        candidate_score = rep(-Inf, length(candidate_indices))

        for (candidate_index in seq_along(candidate_indices)) {
          triangle_indices = candidate_indices[[candidate_index]]
          score_1 = triangle_orientation_score(
            vertices = vertices,
            normals = normals,
            triangle = triangle_indices[1, ]
          )
          score_2 = triangle_orientation_score(
            vertices = vertices,
            normals = normals,
            triangle = triangle_indices[2, ]
          )
          if (is.finite(score_1) && is.finite(score_2)) {
            candidate_score[candidate_index] = min(score_1, score_2) +
              0.01 * (score_1 + score_2)
          }
        }

        best_index = which.max(candidate_score)
        indices[counter:(counter + 1L), ] = candidate_indices[[best_index]]
        norm_indices[counter:(counter + 1L), ] = candidate_indices[[best_index]]
        tex_indices[counter:(counter + 1L), ] = candidate_tex[[best_index]]
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
triangle_orientation_score = function(vertices, normals, triangle) {
  vertex_ids = triangle + 1L
  edge1 = vertices[vertex_ids[2], ] - vertices[vertex_ids[1], ]
  edge2 = vertices[vertex_ids[3], ] - vertices[vertex_ids[1], ]
  face_normal = cross(edge1, edge2)
  if (vector_length(face_normal) <= 1e-8) {
    return(-Inf)
  }

  reference_normal = colMeans(normals[vertex_ids, , drop = FALSE])
  if (vector_length(reference_normal) <= 1e-8) {
    return(-Inf)
  }

  return(sum(face_normal * reference_normal))
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

  counter = 1L
  for (vertex_index in 0:(ring_size - 1L)) {
    next_index = (vertex_index + 1L) %% ring_size

    indices[counter, ] = c(start_center, next_index, vertex_index)
    norm_indices[counter, ] = c(start_normal, start_normal, start_normal)
    tex_indices[counter, ] = c(start_tex, start_tex, start_tex)
    counter = counter + 1L

    indices[counter, ] = c(
      end_center,
      last_ring_vertex + vertex_index,
      last_ring_vertex + next_index
    )
    norm_indices[counter, ] = c(end_normal, end_normal, end_normal)
    tex_indices[counter, ] = c(end_tex, end_tex, end_tex)
    counter = counter + 1L
  }

  return(list(
    indices = indices,
    norm_indices = norm_indices,
    tex_indices = tex_indices
  ))
}

#' @keywords internal
ribbon_texture_u_max = function() {
  1 - 1e-6
}

#' @keywords internal
build_single_ribbon_chain_mesh = function(
  frame,
  ribbon_widths,
  ribbon_thicknesses,
  profile_exponents,
  cross_section_resolution,
  chain_id,
  chain_index
) {
  profile = ribbon_cross_section_profile(
    ribbon_width = ribbon_widths[1],
    ribbon_thickness = ribbon_thicknesses[1],
    cross_section_resolution = cross_section_resolution,
    shape_exponent = profile_exponents[1]
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
    u_values = pmin(u_values, ribbon_texture_u_max())
  }

  for (ring_index in seq_len(ring_count)) {
    profile = ribbon_cross_section_profile(
      ribbon_width = ribbon_widths[ring_index],
      ribbon_thickness = ribbon_thicknesses[ring_index],
      cross_section_resolution = cross_section_resolution,
      shape_exponent = profile_exponents[ring_index]
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
  texcoords[end_center_tex, ] = c(ribbon_texture_u_max(), 0.5)

  side_data = connect_ribbon_rings(
    vertices = vertices[seq_len(ring_count * ring_size), , drop = FALSE],
    ring_count = ring_count,
    ring_size = ring_size,
    normals = normals[seq_len(ring_count * ring_size), , drop = FALSE]
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
  cross_section_resolution = 24,
  shape_exponent = 4
) {
  cross_section_resolution = normalize_cross_section_resolution(
    cross_section_resolution
  )
  if (
    length(shape_exponent) != 1L ||
      !is.finite(shape_exponent) ||
      shape_exponent <= 0
  ) {
    stop("shape_exponent must be a positive number")
  }
  half_width = ribbon_width / 2
  half_thickness = ribbon_thickness / 2
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
  if (
    length(cross_section_resolution) != 1L || is.na(cross_section_resolution)
  ) {
    stop(
      "cross_section_resolution must be an integer greater than or equal to 8"
    )
  }

  cross_section_resolution = as.integer(cross_section_resolution)
  if (
    !is.finite(cross_section_resolution) ||
      cross_section_resolution < 8L
  ) {
    stop(
      "cross_section_resolution must be an integer greater than or equal to 8"
    )
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

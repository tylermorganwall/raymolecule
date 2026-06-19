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
#' @param model_id Default `NA_integer_`. PDB `MODEL` identifier(s) to render.
#'   The default `NA` renders all parsed models in the ensemble.
#' @param ribbon_width Default `2`. Width of the ribbon cross-section.
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
#' @param material Default `rayrender::glossy`. Optional rayrender material
#'   used to initialize the mesh material when `material_vertex` is not
#'   supplied. Must be either `glossy`, `diffuse`, or `dielectric`.
#' @param material_args Default `list()`. Named list of additional arguments
#'   passed to `material`. Arguments supplied by raymolecule for colors and
#'   textures override entries with the same names. For example, use
#'   `list(gloss = 0.35, reflectance = 0.12)` with `rayrender::glossy`, or
#'   `list(sigma = 0.4)` with `rayrender::diffuse`.
#' @param material_vertex Default `rayvertex::material_list(type = "phong")`.
#'   Mesh material template.
#' @param raster_ambient_mix Default `0.5`. Fraction of raster ribbon color
#'   contributed by the ambient, unlit material channel. Values closer to `0`
#'   emphasize diffuse directional lighting; values closer to `1` flatten
#'   lighting and preserve texture/color more directly.
#' @param show_hetero_atoms Default `TRUE`. If `TRUE`, display non-protein
#'   `HETATM` records as bare spheres alongside the ribbon.
#' @param show_hetero_bonds Default `TRUE`. If `TRUE`, display bonds between
#'   shown hetero atoms as a ball-and-stick overlay alongside the ribbon.
#' @param show_waters Default `FALSE`. If `TRUE`, include water `HETATM`
#'   records when `show_hetero_atoms = TRUE` and `show_hetero_bonds = TRUE`.
#' @param show_protein_atoms Default `FALSE`. If `TRUE`, display protein
#'   `ATOM` records as small spheres alongside the ribbon.
#' @param show_protein_bonds Default `FALSE`. If `TRUE`, display inferred
#'   covalent bonds between protein `ATOM` records as a thin stick overlay.
#' @param atom_scale Default `1`. Multiplier applied to the radii of optional
#'   atom overlays.
#' @param bond_width Default `1`. Multiplier applied to the radii of optional
#'   bond overlays.
#' @param use_vertex_normals Default `FALSE`. If `TRUE`, attach the ribbon
#'   mesh's swept vertex normals. If `FALSE`, scenes omit explicit vertex
#'   normals and pathtraced renders use flat face normals.
#' @param verbose Default `FALSE`. If `TRUE`, report the PDB name and model
#'   identifiers being rendered.
#'
#' @return Raymesh scene containing a ribbon mesh.
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' ribbon_file = download_pdb("2w5o", out_dir = tempdir(), overwrite = TRUE)
#' ribbon_model = read_pdb(ribbon_file, verbose = TRUE)
#'
#' # Start with a centered raster ribbon using the default ribbon width,
#' # thickness, color mode, and ligand overlays.
#' ribbon_model |>
#'   generate_ribbon_scene() |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey12"
#'   )
#'
#' # This pathtraced version widens the ribbon, increases cross-section
#' # resolution, applies a custom UV texture, turns on atom/bond overlays, and
#' # uses the mesh vertex normals.
#' texture_file = tempfile(fileext = ".png")
#' grDevices::png(texture_file, width = 64, height = 8, bg = "transparent")
#' graphics::par(mar = c(0, 0, 0, 0))
#' graphics::image(
#'   matrix(seq(0, 1, length.out = 64), ncol = 1),
#'   col = grDevices::hcl.colors(64, "Spectral", rev = TRUE),
#'   axes = FALSE,
#'   xlab = "",
#'   ylab = ""
#' )
#' grDevices::dev.off()
#'
#' ribbon_model |>
#'   generate_ribbon_scene(
#'     x = 0,
#'     y = 0,
#'     z = 0,
#'     scale = 1,
#'     center = TRUE,
#'     ribbon_width = 1.8,
#'     ribbon_thickness = 0.3,
#'     cross_section_resolution = 32,
#'     subdivisions = 10,
#'     color_mode = "uv",
#'     texture = texture_file,
#'     material = rayrender::diffuse,
#'     material_args = list(sigma = 0.2),
#'     show_hetero_atoms = TRUE,
#'     show_hetero_bonds = TRUE,
#'     show_waters = TRUE,
#'     show_protein_atoms = TRUE,
#'     show_protein_bonds = TRUE,
#'     atom_scale = 1.2,
#'     bond_width = 0.8,
#'     use_vertex_normals = TRUE,
#'     verbose = TRUE
#'   ) |>
#'   render_model(pathtrace = TRUE, width = 800, height = 800, samples = 32)
#'
#' # Start the beta-barrel example with the default UV texture in raster mode.
#' # We rotate the model to match the render in the protein data bank.
#' barrel_model = read_pdb(download_pdb("4fsp", out_dir = tempdir()))
#' barrel_scene = barrel_model |>
#'   generate_ribbon_scene(color_mode = "uv", raster_ambient_mix = 0.7)
#' barrel_render = render_model(
#'   barrel_scene,
#'   pathtrace = FALSE,
#'   width = 500,
#'   height = 500,
#'   background = "white",
#'   lookfrom = c(-89.95, 66.11, -109.95),
#'   angle = c(-60, 270, 180),
#'   lookat = c(6.06, -7.62, 1.41),
#'   fov = 27.2,
#'   plot = FALSE
#' )
#'
#' pdb_4fsp_image = rayimage::render_title(system.file(
#'   "extdata",
#'   "4fsp_assembly-1.jpeg",
#'   package = "raymolecule",
#'   mustWork = TRUE
#' ), title_text = "Protein Data Bank Render")
#'
#' ray_4fsp_image = rayimage::render_title(
#'   barrel_render,
#'   title_text = "Raymolecule Render"
#' )
#' rayimage::plot_image_grid(list(pdb_4fsp_image, ray_4fsp_image),dim=c(1,2))
#' # Raising the raster ambient mix gives the same barrel less directional
#' # lighting
#' barrel_model |>
#'   generate_ribbon_scene(
#'     color_mode = "uv",
#'     texture = NULL,
#'     raster_ambient_mix = 0.8
#'   ) |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey12"
#'   )
#'
#' # Start the multi-chain example with the automatic chain color palette.
#' multi_chain_model = read_pdb(download_pdb("1xn1", out_dir = tempdir()))
#' multi_chain_model |>
#'   generate_ribbon_scene(color_mode = "chain") |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey80"
#'   )
#'
#' # An explicit chain color map replaces the automatic palette without
#' # changing the ribbon geometry.
#' chain_ids = unique(multi_chain_model$residues$chain_id)
#' chain_ids = chain_ids[!is.na(chain_ids)]
#' chain_colors = stats::setNames(grDevices::rainbow(length(chain_ids)), chain_ids)
#' multi_chain_model |>
#'   generate_ribbon_scene(
#'     color_mode = "chain",
#'     chain_colors = chain_colors,
#'     material_vertex = rayvertex::material_list(type = "phong")
#'   ) |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey80"
#'   )
#'
#' # Start the NMR ensemble view by rendering every parsed model.
#' ensemble_model = read_pdb(download_pdb("1co1", out_dir = tempdir()))
#' ensemble_model |>
#'   generate_ribbon_scene(center = FALSE) |>
#'   render_model(
#'     pathtrace = TRUE,
#'     fov=26,
#'		 lookfrom = c(-45.95, 58.56, 79.95),
#'		 lookat = c(-7.19, 3.87, -1.52) ,
#'     width = 800,
#'     height = 800,
#'     background = "black"
#'   )
#'
#' # Selecting three model IDs shows a cleaner subset of the same ensemble.
#' ensemble_model |>
#'	 generate_ribbon_scene(model_id = c(1L, 5L, 10L), center = FALSE) |>
#'	 render_model(
#'		 pathtrace = TRUE,
#'     fov=26,
#'		 lookfrom = c(-45.95, 58.56, 79.95),
#'		 lookat = c(-7.19, 3.87, -1.52) ,
#'		 width = 800,
#'		 height = 800,
#'		 background = "black"
#'	 )
generate_ribbon_scene = function(
  model,
  x = 0,
  y = 0,
  z = 0,
  scale = 1,
  center = TRUE,
  model_id = NA_integer_,
  ribbon_width = 2,
  ribbon_thickness = 0.25,
  cross_section_resolution = 24,
  subdivisions = 8,
  color_mode = c("chain", "uv"),
  chain_colors = NULL,
  texture = NULL,
  material = rayrender::glossy,
  material_args = list(),
  material_vertex = rayvertex::material_list(type = "phong"),
  raster_ambient_mix = 0.5,
  show_hetero_atoms = TRUE,
  show_hetero_bonds = TRUE,
  show_waters = FALSE,
  show_protein_atoms = FALSE,
  show_protein_bonds = FALSE,
  atom_scale = 1,
  bond_width = 1,
  use_vertex_normals = FALSE,
  verbose = FALSE
) {
  if (!identical(model$pdb_type, "pdb") || !is.data.frame(model$residues)) {
    stop("generate_ribbon_scene() requires a PDB model produced by read_pdb()")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be a positive number")
  }
  if (length(atom_scale) != 1L || !is.finite(atom_scale) || atom_scale <= 0) {
    stop("atom_scale must be a positive number")
  }
  if (length(bond_width) != 1L || !is.finite(bond_width) || bond_width <= 0) {
    stop("bond_width must be a positive number")
  }
  if (
    length(raster_ambient_mix) != 1L ||
      !is.finite(raster_ambient_mix) ||
      raster_ambient_mix < 0 ||
      raster_ambient_mix > 1
  ) {
    stop("raster_ambient_mix must be a finite number between 0 and 1")
  }
  model = select_pdb_models(model, model_id)
  if (verbose) {
    message(sprintf(
      "Rendering %s PDB models %s",
      pdb_display_name(model),
      format_pdb_model_set(pdb_model_labels(model))
    ))
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
  if (
    !is.logical(use_vertex_normals) ||
      length(use_vertex_normals) != 1L ||
      is.na(use_vertex_normals)
  ) {
    stop("use_vertex_normals must be TRUE or FALSE")
  }
  if (
    !is.logical(show_hetero_atoms) ||
      length(show_hetero_atoms) != 1L ||
      is.na(show_hetero_atoms)
  ) {
    stop("show_hetero_atoms must be TRUE or FALSE")
  }
  if (
    !is.logical(show_hetero_bonds) ||
      length(show_hetero_bonds) != 1L ||
      is.na(show_hetero_bonds)
  ) {
    stop("show_hetero_bonds must be TRUE or FALSE")
  }
  if (
    !is.logical(show_waters) || length(show_waters) != 1L || is.na(show_waters)
  ) {
    stop("show_waters must be TRUE or FALSE")
  }
  if (
    !is.logical(show_protein_atoms) ||
      length(show_protein_atoms) != 1L ||
      is.na(show_protein_atoms)
  ) {
    stop("show_protein_atoms must be TRUE or FALSE")
  }
  if (
    !is.logical(show_protein_bonds) ||
      length(show_protein_bonds) != 1L ||
      is.na(show_protein_bonds)
  ) {
    stop("show_protein_bonds must be TRUE or FALSE")
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
  visible_hetero_atoms = select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = show_hetero_atoms,
    show_waters = show_waters,
    show_protein_atoms = FALSE
  )
  visible_protein_atoms = select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = FALSE,
    show_waters = FALSE,
    show_protein_atoms = show_protein_atoms
  )
  visible_protein_bond_atoms = select_ribbon_display_atoms(
    model = model,
    show_hetero_atoms = FALSE,
    show_waters = FALSE,
    show_protein_atoms = show_protein_bonds
  )
  visible_hetero_bonds = select_ribbon_display_bonds(
    model = model,
    atoms = visible_hetero_atoms,
    show_hetero_bonds = show_hetero_bonds,
    show_protein_bonds = FALSE
  )
  visible_protein_bonds = select_ribbon_display_bonds(
    model = model,
    atoms = visible_protein_bond_atoms,
    show_hetero_bonds = FALSE,
    show_protein_bonds = show_protein_bonds
  )

  use_material = !missing(material) && missing(material_vertex)
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
      material = material,
      material_vertex = material_vertex,
      use_material = use_material,
      material_args = material_args,
      color_mode = color_mode,
      chain_id = chain_mesh$chain_id,
      chain_color = if (is.null(chain_color_lookup)) {
        NULL
      } else {
        chain_color_lookup[chain_mesh$chain_id][[1]]
      },
      texture = texture,
      raster_ambient_mix = raster_ambient_mix
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
  if (nrow(visible_hetero_bonds) > 0L) {
    mesh_list = c(
      mesh_list,
      build_ribbon_bond_meshes(
        atoms = visible_hetero_atoms,
        bonds = visible_hetero_bonds,
        center_shift = center_shift,
        scale = scale,
        offset = offset,
        material_vertex = material_vertex,
        bond_radius_scale = 1 * bond_width
      )
    )
  }
  if (nrow(visible_protein_bonds) > 0L) {
    mesh_list = c(
      mesh_list,
      build_ribbon_bond_meshes(
        atoms = visible_protein_bond_atoms,
        bonds = visible_protein_bonds,
        center_shift = center_shift,
        scale = scale,
        offset = offset,
        material_vertex = material_vertex,
        bond_radius_scale = 0.2 * bond_width
      )
    )
  }
  if (nrow(visible_hetero_atoms) > 0L) {
    mesh_list = c(
      mesh_list,
      build_ribbon_atom_meshes(
        atoms = visible_hetero_atoms,
        center_shift = center_shift,
        scale = scale,
        offset = offset,
        material_vertex = material_vertex,
        atom_radius_scale = 0.45 * atom_scale
      )
    )
  }
  if (nrow(visible_protein_atoms) > 0L) {
    mesh_list = c(
      mesh_list,
      build_ribbon_atom_meshes(
        atoms = visible_protein_atoms,
        center_shift = center_shift,
        scale = scale,
        offset = offset,
        material_vertex = material_vertex,
        atom_radius_scale = 0.18 * atom_scale
      )
    )
  }

  mesh_scene = rayvertex::scene_from_list(mesh_list)
  return(mesh_scene)
}

#' @keywords internal
select_pdb_models = function(model, model_id = NA_integer_) {
  if (length(model_id) == 1L && is.na(model_id)) {
    return(model)
  }
  if (
    length(model_id) == 0L ||
      any(is.na(model_id)) ||
      any(!is.finite(model_id))
  ) {
    stop("model_id must be NA or a finite PDB MODEL identifier")
  }
  model_id = as.integer(model_id)

  available_model_ids = pdb_model_labels(model)
  available_model_ids = available_model_ids[!is.na(available_model_ids)]
  missing_model_ids = setdiff(model_id, available_model_ids)
  if (length(missing_model_ids) > 0L) {
    stop(sprintf(
      "PDB model_id value(s) not found: %s. Available PDB models: %s",
      paste(missing_model_ids, collapse = ", "),
      format_pdb_model_set(available_model_ids)
    ))
  }

  for (name in names(model)) {
    value = model[[name]]
    if (!is.data.frame(value) || !"model_index" %in% names(value)) {
      next
    }
    model[[name]] = value[
      pdb_model_row_labels(value) %in% model_id,
      ,
      drop = FALSE
    ]
    rownames(model[[name]]) = NULL
  }

  model
}

#' @keywords internal
pdb_model_row_labels = function(data) {
  if (!"model_index" %in% names(data)) {
    return(rep(NA_integer_, nrow(data)))
  }
  if (!"model_id" %in% names(data) || all(is.na(data$model_id))) {
    return(data$model_index)
  }
  ifelse(is.na(data$model_id), data$model_index, data$model_id)
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
    if (
      !is.matrix(fallback_norm_indices) || ncol(fallback_norm_indices) != 3L
    ) {
      stop(
        "build_flat_mesh_normals() fallback_norm_indices must be an m x 3 matrix"
      )
    }
    if (nrow(fallback_norm_indices) != nrow(indices)) {
      stop("build_flat_mesh_normals() fallback_norm_indices must match indices")
    }
  }

  triangle_count = nrow(indices)
  normals = matrix(NA_real_, nrow = triangle_count, ncol = 3L)
  norm_indices = matrix(
    rep(seq_len(triangle_count) - 1L, each = 3L),
    ncol = 3L,
    byrow = TRUE
  )
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
      if (
        !is.null(previous_normal) && sum(fallback_normal * previous_normal) < 0
      ) {
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
  show_waters,
  show_protein_atoms = FALSE
) {
  if (
    (!show_hetero_atoms && !show_protein_atoms) ||
      is.null(model$atoms) ||
      nrow(model$atoms) == 0L
  ) {
    return(empty_pdb_atoms())
  }

  atom_rows = list()
  if (show_hetero_atoms) {
    hetero_atoms = model$atoms[model$atoms$record == "HETATM", , drop = FALSE]
    if (!show_waters) {
      hetero_atoms = hetero_atoms[
        !is_water_residue_name(hetero_atoms$res_name),
        ,
        drop = FALSE
      ]
    }
    atom_rows[[length(atom_rows) + 1L]] = hetero_atoms
  }
  if (show_protein_atoms) {
    atom_rows[[length(atom_rows) + 1L]] = model$atoms[
      model$atoms$record == "ATOM",
      ,
      drop = FALSE
    ]
  }
  return(bind_pdb_rows(atom_rows, empty_pdb_atoms()))
}

#' @keywords internal
select_ribbon_display_bonds = function(
  model,
  atoms,
  show_hetero_bonds,
  show_protein_bonds = FALSE
) {
  if (
    (!show_hetero_bonds && !show_protein_bonds) ||
      nrow(atoms) == 0L
  ) {
    return(empty_pdb_bonds())
  }

  bond_rows = list()
  if (show_hetero_bonds && !is.null(model$bonds) && nrow(model$bonds) > 0L) {
    visible_hetero_ids = atoms$index[atoms$record == "HETATM"]
    hetero_bonds = model$bonds[
      model$bonds$from %in%
        visible_hetero_ids &
        model$bonds$to %in% visible_hetero_ids,
      ,
      drop = FALSE
    ]
    if (nrow(hetero_bonds) > 0L) {
      bond_rows[[length(bond_rows) + 1L]] = hetero_bonds
    }
  }

  if (show_protein_bonds) {
    inferred_bonds = infer_pdb_atom_bonds(
      atoms[atoms$record == "ATOM", , drop = FALSE]
    )
    if (nrow(inferred_bonds) > 0L) {
      bond_rows[[length(bond_rows) + 1L]] = inferred_bonds
    }
  }

  bonds = bind_pdb_rows(bond_rows, empty_pdb_bonds())
  return(normalize_ribbon_display_bonds(bonds))
}

#' @keywords internal
infer_pdb_atom_bonds = function(atoms) {
  if (!is.data.frame(atoms) || nrow(atoms) < 2L) {
    return(empty_pdb_bonds())
  }

  atom_order = if ("atom_order" %in% names(atoms)) {
    atoms$atom_order
  } else {
    atoms$index
  }
  atoms = atoms[
    order(
      atoms$model_index,
      atoms$chain_id,
      atoms$res_seq,
      atoms$i_code,
      atom_order
    ),
  ]
  bond_rows = list()
  residue_keys = paste(
    atoms$model_index,
    atoms$chain_id,
    atoms$res_seq,
    atoms$i_code,
    sep = "\r"
  )
  residue_key_order = unique(residue_keys)

  for (residue_key in residue_key_order) {
    residue_atoms = atoms[residue_keys == residue_key, , drop = FALSE]
    if (nrow(residue_atoms) < 2L) {
      next
    }
    for (i in seq_len(nrow(residue_atoms) - 1L)) {
      for (j in (i + 1L):nrow(residue_atoms)) {
        if (
          !pdb_atom_alt_locs_compatible(
            residue_atoms$alt_loc[i],
            residue_atoms$alt_loc[j]
          )
        ) {
          next
        }
        if (identical(residue_atoms$atom_name[i], residue_atoms$atom_name[j])) {
          next
        }
        if (pdb_atom_pair_is_bonded(residue_atoms[i, ], residue_atoms[j, ])) {
          bond_rows[[length(bond_rows) + 1L]] = pdb_inferred_bond_row(
            residue_atoms[i, ],
            residue_atoms[j, ]
          )
        }
      }
    }
  }

  peptide_bonds = infer_pdb_peptide_bonds(atoms, residue_keys)
  if (nrow(peptide_bonds) > 0L) {
    bond_rows[[length(bond_rows) + 1L]] = peptide_bonds
  }

  disulfide_bonds = infer_pdb_disulfide_bonds(atoms)
  if (nrow(disulfide_bonds) > 0L) {
    bond_rows[[length(bond_rows) + 1L]] = disulfide_bonds
  }

  bonds = bind_pdb_rows(bond_rows, empty_pdb_bonds())
  normalize_ribbon_display_bonds(bonds)
}

#' @keywords internal
infer_pdb_peptide_bonds = function(atoms, residue_keys) {
  residue_key_order = unique(residue_keys)
  if (length(residue_key_order) < 2L) {
    return(empty_pdb_bonds())
  }

  bond_rows = list()
  residue_starts = match(residue_key_order, residue_keys)
  residues = atoms[residue_starts, , drop = FALSE]

  for (i in seq_len(nrow(residues) - 1L)) {
    if (
      residues$model_index[i] != residues$model_index[i + 1L] ||
        !identical(residues$chain_id[i], residues$chain_id[i + 1L])
    ) {
      next
    }

    left_atoms = atoms[residue_keys == residue_key_order[i], , drop = FALSE]
    right_atoms = atoms[
      residue_keys == residue_key_order[i + 1L],
      ,
      drop = FALSE
    ]
    c_atoms = left_atoms[left_atoms$atom_name == "C", , drop = FALSE]
    n_atoms = right_atoms[right_atoms$atom_name == "N", , drop = FALSE]
    if (nrow(c_atoms) == 0L || nrow(n_atoms) == 0L) {
      next
    }

    for (c_index in seq_len(nrow(c_atoms))) {
      for (n_index in seq_len(nrow(n_atoms))) {
        if (
          !pdb_atom_alt_locs_compatible(
            c_atoms$alt_loc[c_index],
            n_atoms$alt_loc[n_index]
          )
        ) {
          next
        }
        if (pdb_atom_pair_is_bonded(c_atoms[c_index, ], n_atoms[n_index, ])) {
          bond_rows[[length(bond_rows) + 1L]] = pdb_inferred_bond_row(
            c_atoms[c_index, ],
            n_atoms[n_index, ]
          )
        }
      }
    }
  }

  bind_pdb_rows(bond_rows, empty_pdb_bonds())
}

#' @keywords internal
infer_pdb_disulfide_bonds = function(atoms) {
  sulfur_atoms = atoms[
    atoms$atom_name == "SG" &
      toupper(atoms$res_name) == "CYS",
    ,
    drop = FALSE
  ]
  if (nrow(sulfur_atoms) < 2L) {
    return(empty_pdb_bonds())
  }

  bond_rows = list()
  for (i in seq_len(nrow(sulfur_atoms) - 1L)) {
    for (j in (i + 1L):nrow(sulfur_atoms)) {
      if (sulfur_atoms$model_index[i] != sulfur_atoms$model_index[j]) {
        next
      }
      if (
        identical(sulfur_atoms$chain_id[i], sulfur_atoms$chain_id[j]) &&
          sulfur_atoms$res_seq[i] == sulfur_atoms$res_seq[j] &&
          identical(sulfur_atoms$i_code[i], sulfur_atoms$i_code[j])
      ) {
        next
      }
      if (
        !pdb_atom_alt_locs_compatible(
          sulfur_atoms$alt_loc[i],
          sulfur_atoms$alt_loc[j]
        )
      ) {
        next
      }
      if (pdb_atom_pair_is_bonded(sulfur_atoms[i, ], sulfur_atoms[j, ])) {
        bond_rows[[length(bond_rows) + 1L]] = pdb_inferred_bond_row(
          sulfur_atoms[i, ],
          sulfur_atoms[j, ]
        )
      }
    }
  }

  bind_pdb_rows(bond_rows, empty_pdb_bonds())
}

#' @keywords internal
pdb_atom_pair_is_bonded = function(atom1, atom2) {
  distance = sqrt(
    (atom1$x - atom2$x)^2 +
      (atom1$y - atom2$y)^2 +
      (atom1$z - atom2$z)^2
  )
  if (!is.finite(distance) || distance <= 0.4) {
    return(FALSE)
  }

  max_distance =
    (pdb_covalent_radius(atom1$type) + pdb_covalent_radius(atom2$type)) *
    1.25 +
    0.15
  distance <= max_distance
}

#' @keywords internal
pdb_covalent_radius = function(element) {
  element = normalize_atom_element_symbol(element)
  radius = suppressWarnings(tryCatch(
    PeriodicTable::rcov(element),
    error = function(e) NA_real_
  ))
  if (length(radius) != 1L || is.na(radius) || !is.finite(radius)) {
    return(0.77)
  }
  radius
}

#' @keywords internal
pdb_atom_alt_locs_compatible = function(alt_loc1, alt_loc2) {
  alt_loc1 = trimws(alt_loc1)
  alt_loc2 = trimws(alt_loc2)
  if (is.na(alt_loc1)) {
    alt_loc1 = ""
  }
  if (is.na(alt_loc2)) {
    alt_loc2 = ""
  }
  !nzchar(alt_loc1) || !nzchar(alt_loc2) || identical(alt_loc1, alt_loc2)
}

#' @keywords internal
pdb_inferred_bond_row = function(atom1, atom2) {
  data.frame(
    from = atom1$index,
    to = atom2$index,
    number = 1L,
    model_id = atom1$model_id,
    model_index = atom1$model_index,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
normalize_ribbon_display_bonds = function(bonds) {
  if (nrow(bonds) == 0L) {
    return(empty_pdb_bonds())
  }

  if (!"model_id" %in% names(bonds)) {
    bonds$model_id = NA_integer_
  }
  if (!"model_index" %in% names(bonds)) {
    bonds$model_index = NA_integer_
  }
  bonds = data.frame(
    from = pmin(bonds$from, bonds$to),
    to = pmax(bonds$from, bonds$to),
    number = bonds$number,
    model_id = bonds$model_id,
    model_index = bonds$model_index,
    stringsAsFactors = FALSE
  )
  bonds = bonds[bonds$from != bonds$to, , drop = FALSE]
  if (nrow(bonds) == 0L) {
    return(empty_pdb_bonds())
  }

  bond_key = paste(
    bonds$from,
    bonds$to,
    ifelse(is.na(bonds$model_index), "", bonds$model_index),
    sep = "\r"
  )
  max_number = ave(bonds$number, bond_key, FUN = max)
  keep_bond = !duplicated(bond_key)
  bonds = bonds[keep_bond, , drop = FALSE]
  bonds$number = max_number[keep_bond]
  bonds = bonds[order(bonds$model_index, bonds$from, bonds$to), , drop = FALSE]
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
  return(
    min(
      ribbon_display_atom_radius(element1),
      ribbon_display_atom_radius(element2)
    ) *
      scale
  )
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
  material_vertex,
  bond_radius_scale = 0.6
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
    radius = ribbon_display_bond_radius(
      from_atom$type,
      to_atom$type,
      scale = bond_radius_scale
    )

    from_material = material_vertex
    from_color = display_atom_color(from_atom$type)
    from_material$diffuse = convert_color(from_color)
    from_material$ambient = convert_color(from_color)
    from_material$ambient_intensity = max(0.1, from_material$ambient_intensity)
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
  material_vertex,
  atom_radius_scale = 0.45
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
      radius = ribbon_display_atom_radius(
        atoms$type[i],
        scale = atom_radius_scale
      ),
      low_poly = FALSE,
      material = atom_material
    )
  }

  return(meshes)
}

#' @keywords internal
prepare_ribbon_material = function(
  material,
  material_vertex,
  use_material,
  material_args,
  color_mode,
  chain_id,
  chain_color,
  texture,
  raster_ambient_mix
) {
  if (use_material) {
    package_args = list()
    if (identical(color_mode, "chain") && !is.null(chain_color)) {
      package_args$color = chain_color
      package_args$image_texture = ""
    } else if (identical(color_mode, "uv") && !is.null(texture)) {
      package_args$color = "white"
      package_args$image_texture = texture
    }
    mesh_material = rayrender_material_to_vertex_material(
      material,
      material_args = material_args,
      package_args = package_args
    )
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
    mesh_material$diffuse_texname = texture
    mesh_material$ambient_texname = texture
  }

  mesh_material$ambient_intensity = raster_ambient_mix
  mesh_material$diffuse_intensity = 1 - raster_ambient_mix

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
      graphics::plot.window(
        xlim = c(0, 1),
        ylim = c(0, 1),
        xaxs = "i",
        yaxs = "i"
      )
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
resolve_mesh_material = function(
  material,
  material_vertex,
  use_material,
  material_args = list(),
  package_args = list()
) {
  if (isTRUE(use_material)) {
    return(rayrender_material_to_vertex_material(
      material,
      material_args = material_args,
      package_args = package_args
    ))
  }
  return(material_vertex)
}

#' @keywords internal
rayrender_material_to_vertex_material = function(
  material,
  material_args = list(),
  package_args = list()
) {
  if (!is.function(material)) {
    stop("material must be a rayrender material function")
  }

  resolved_material_args = resolve_rayrender_material_args(
    material = material,
    material_args = material_args,
    package_args = package_args
  )
  material_info = do.call(material, resolved_material_args)
  if (!inherits(material_info, "ray_material") || length(material_info) == 0L) {
    stop("material() must return a rayrender material")
  }

  material_info = material_info[[1]]
  if (!material_info$type %in% c(7L, 4L, 1L, 3L)) {
    stop("material() must be either `glossy`, `diffuse`, or `dielectric`")
  }

  base_color = material_info$properties[[1]][1:3]
  mesh_material = rayvertex::material_list(
    diffuse = base_color,
    ambient = base_color,
    ambient_intensity = 0.2,
    type = "phong"
  )

  if (material_info$type %in% c(4L, 1L)) {
    mesh_material$type = "diffuse"
    mesh_material$specular = c(0, 0, 0)
    mesh_material$specular_intensity = 0
    if (length(material_info$sigma) == 1L && is.finite(material_info$sigma)) {
      mesh_material$sigma = material_info$sigma
    }
  } else if (material_info$type == 7L) {
    glossy_info = material_info$glossyinfo[[1]]
    if (length(glossy_info) >= 6L && all(is.finite(glossy_info[4:6]))) {
      reflectance = mean(glossy_info[4:6])
      mesh_material$reflection_intensity = reflectance
      mesh_material$specular_intensity = reflectance
    }
    if (length(glossy_info) >= 3L && all(is.finite(glossy_info[2:3]))) {
      gloss = 1 - 2 * sqrt(mean(glossy_info[2:3]))
      gloss = max(0, min(1, gloss))
      mesh_material$reflection_sharpness = gloss
      mesh_material$shininess = max(1, 100 * gloss)
    }
  } else if (material_info$type == 3L) {
    mesh_material$type = "phong"
    if (length(material_info$properties[[1]]) >= 4L) {
      mesh_material$ior = material_info$properties[[1]][4]
    }
  }

  if (!is.null(material_info$image) && nzchar(material_info$image)) {
    mesh_material$diffuse_texname = material_info$image
  }
  if (
    !is.null(material_info$bump_texture) && nzchar(material_info$bump_texture)
  ) {
    mesh_material$bump_texname = material_info$bump_texture
  }
  if (
    !is.null(material_info$bump_intensity) &&
      length(material_info$bump_intensity) == 1L &&
      is.finite(material_info$bump_intensity)
  ) {
    mesh_material$bump_intensity = material_info$bump_intensity
  }

  mesh_material = attach_rayrender_material(
    mesh_material = mesh_material,
    material = material,
    material_args = material_args,
    package_args = package_args
  )
  return(mesh_material)
}

#' @keywords internal
attach_rayrender_material = function(
  mesh_material,
  material,
  material_args = list(),
  package_args = list()
) {
  material_args = validate_rayrender_material_args(
    material_args,
    "material_args"
  )
  package_args = validate_rayrender_material_args(package_args, "package_args")
  attr(mesh_material, "raymolecule_rayrender_material") = do.call(
    material,
    resolve_rayrender_material_args(
      material = material,
      material_args = material_args,
      package_args = package_args
    )
  )
  attr(mesh_material, "raymolecule_rayrender_material_function") = material
  attr(mesh_material, "raymolecule_rayrender_material_args") = material_args
  attr(mesh_material, "raymolecule_rayrender_package_args") = package_args
  return(mesh_material)
}

#' @keywords internal
update_rayrender_material_package_args = function(
  mesh_material,
  package_args = list()
) {
  material = attr(mesh_material, "raymolecule_rayrender_material_function")
  material_args = attr(mesh_material, "raymolecule_rayrender_material_args")
  if (is.null(material) || is.null(material_args)) {
    return(mesh_material)
  }

  attach_rayrender_material(
    mesh_material = mesh_material,
    material = material,
    material_args = material_args,
    package_args = package_args
  )
}

#' @keywords internal
resolve_rayrender_material_args = function(
  material,
  material_args = list(),
  package_args = list()
) {
  material_args = validate_rayrender_material_args(
    material_args,
    "material_args"
  )
  package_args = validate_rayrender_material_args(package_args, "package_args")

  material_formals = names(formals(material))
  if (is.null(material_formals)) {
    material_formals = character()
  }
  accepts_dots = "..." %in% material_formals

  if (!accepts_dots) {
    unknown_args = setdiff(names(material_args), material_formals)
    if (length(unknown_args) > 0L) {
      stop(sprintf(
        "Unknown rayrender material argument(s): %s",
        paste(unknown_args, collapse = ", ")
      ))
    }
    package_args = package_args[names(package_args) %in% material_formals]
  }

  if (length(package_args) > 0L) {
    material_args[names(package_args)] = package_args
  }
  return(material_args)
}

#' @keywords internal
validate_rayrender_material_args = function(material_args, argument_name) {
  if (is.null(material_args)) {
    return(list())
  }
  if (!is.list(material_args)) {
    stop(sprintf("%s must be a named list", argument_name))
  }
  if (length(material_args) == 0L) {
    return(material_args)
  }

  arg_names = names(material_args)
  if (is.null(arg_names) || any(is.na(arg_names) | !nzchar(arg_names))) {
    stop(sprintf("%s must be a named list", argument_name))
  }
  if (anyDuplicated(arg_names)) {
    stop(sprintf("%s must not contain duplicate names", argument_name))
  }
  return(material_args)
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

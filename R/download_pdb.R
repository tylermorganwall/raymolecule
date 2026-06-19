#' Download a PDB Structure from RCSB
#'
#' Looks up a protein name in the RCSB PDB search API and downloads the best
#' matching legacy PDB file. A 4-character PDB ID can also be supplied directly.
#'
#' @param protein Protein name or 4-character PDB ID.
#' @param out_dir Default `"."`. Directory where the PDB file will be written.
#' @param filename Default `NULL`. Optional output file name. If `NULL`, the
#'   downloaded file is named with the matched PDB ID.
#' @param overwrite Default `FALSE`. Whether to overwrite an existing file.
#' @param max_results Default `10L`. Maximum number of RCSB search matches to
#'   try when downloading a legacy PDB file.
#' @param verbose Default `FALSE`. If `TRUE`, report the matched PDB ID and
#'   destination path.
#'
#' @return Path to the downloaded PDB file.
#' @export
#'
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Start with a direct PDB ID download. A temporary directory, custom
#' # filename, overwrite flag, and verbose output make the file handling clear.
#' pdb_file = download_pdb(
#'   "4fsp",
#'   out_dir = tempdir(),
#'   filename = "outer-membrane-barrel.pdb",
#'   overwrite = TRUE,
#'   verbose = TRUE
#' )
#'
#' read_pdb(pdb_file, verbose = TRUE) |>
#'   generate_ribbon_scene() |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey12"
#'   )
#'
#' # Protein-name lookup searches RCSB and tries up to three legacy PDB matches.
#' hemoglobin_file = download_pdb(
#'   "hemoglobin",
#'   out_dir = tempdir(),
#'   max_results = 3L
#' )
#' read_pdb(hemoglobin_file, verbose = TRUE) |>
#'   generate_ribbon_scene() |>
#'   render_model(
#'     pathtrace = FALSE,
#'     width = 800,
#'     height = 800,
#'     background = "grey12"
#'   )
download_pdb = function(
  protein,
  out_dir = ".",
  filename = NULL,
  overwrite = FALSE,
  max_results = 10L,
  verbose = FALSE
) {
  protein = validate_rcsb_protein_query(protein)
  max_results = validate_rcsb_max_results(max_results)
  verbose = validate_rcsb_flag(verbose, "verbose")
  overwrite = validate_rcsb_flag(overwrite, "overwrite")

  if (is_rcsb_pdb_id(protein)) {
    pdb_ids = toupper(protein)
  } else {
    pdb_ids = rcsb_search_pdb_ids(protein, max_results = max_results)
  }

  if (length(pdb_ids) == 0L) {
    stop(
      sprintf("Cannot find protein '%s' in RCSB PDB", protein),
      call. = FALSE
    )
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(out_dir)) {
    stop(sprintf("Cannot create output directory: %s", out_dir), call. = FALSE)
  }

  failed = character()
  for (pdb_id in pdb_ids) {
    dest = rcsb_pdb_destination(pdb_id, out_dir, filename)
    ok = tryCatch(
      download_rcsb_pdb_file(
        pdb_id = pdb_id,
        dest = dest,
        overwrite = overwrite
      ),
      error = function(e) {
        failed <<- c(
          failed,
          sprintf("%s (%s)", pdb_id, conditionMessage(e))
        )
        FALSE
      }
    )

    if (ok) {
      if (verbose) {
        message(sprintf(
          "Downloaded RCSB PDB %s for '%s' to %s",
          pdb_id,
          protein,
          dest
        ))
      }
      attr(dest, "pdb_id") = pdb_id
      attr(dest, "query") = protein
      attr(dest, "matched_ids") = pdb_ids
      return(dest)
    }
    failed = c(failed, pdb_id)
  }

  stop(
    sprintf(
      paste(
        "RCSB matches were found for protein '%s',",
        "but no legacy PDB file could be downloaded: %s"
      ),
      protein,
      paste(unique(failed), collapse = ", ")
    ),
    call. = FALSE
  )
}

#' @keywords internal
rcsb_search_pdb_ids = function(protein, max_results = 10L) {
  query = build_rcsb_protein_search_query(protein, max_results = max_results)
  response = httr::POST(
    "https://search.rcsb.org/rcsbsearch/v2/query",
    body = query,
    encode = "json",
    httr::accept_json()
  )
  status = httr::status_code(response)

  if (identical(status, 204L)) {
    return(character())
  }
  if (status < 200L || status >= 300L) {
    response_text = httr::content(
      response,
      as = "text",
      encoding = "UTF-8"
    )
    stop(
      sprintf(
        "RCSB search failed with HTTP status %d: %s",
        status,
        response_text
      ),
      call. = FALSE
    )
  }

  content = httr::content(
    response,
    as = "parsed",
    type = "application/json"
  )
  extract_rcsb_pdb_ids(content)
}

#' @keywords internal
build_rcsb_protein_search_query = function(protein, max_results = 10L) {
  list(
    query = list(
      type = "group",
      logical_operator = "and",
      nodes = list(
        list(
          type = "terminal",
          service = "full_text",
          parameters = list(value = protein)
        ),
        list(
          type = "terminal",
          service = "text",
          parameters = list(
            attribute = "rcsb_entry_info.polymer_entity_count_protein",
            operator = "greater",
            value = 0
          )
        )
      )
    ),
    return_type = "entry",
    request_options = list(
      results_content_type = list("experimental"),
      paginate = list(start = 0L, rows = as.integer(max_results)),
      sort = list(list(sort_by = "score", direction = "desc"))
    )
  )
}

#' @keywords internal
extract_rcsb_pdb_ids = function(content) {
  result_set = content$result_set
  if (is.null(result_set)) {
    return(character())
  }
  if (is.character(result_set)) {
    return(unique(toupper(result_set[nzchar(result_set)])))
  }
  if (is.data.frame(result_set) && "identifier" %in% names(result_set)) {
    ids = result_set$identifier
    return(unique(toupper(ids[!is.na(ids) & nzchar(ids)])))
  }

  ids = vapply(
    result_set,
    function(result) {
      if (is.character(result)) {
        return(result[[1]])
      }
      if (!is.null(result$identifier)) {
        return(result$identifier)
      }
      NA_character_
    },
    character(1)
  )
  unique(toupper(ids[!is.na(ids) & nzchar(ids)]))
}

#' @keywords internal
download_rcsb_pdb_file = function(pdb_id, dest, overwrite = FALSE) {
  if (file.exists(dest) && !overwrite) {
    if (is_valid_pdb_download(dest)) {
      return(TRUE)
    }
    stop(sprintf("Destination exists and is not a valid PDB file: %s", dest))
  }

  url = sprintf("https://files.rcsb.org/download/%s.pdb", toupper(pdb_id))
  ok = tryCatch(
    {
      utils::download.file(url, dest, mode = "wb", quiet = TRUE)
      TRUE
    },
    error = function(e) FALSE
  )

  if (!ok || !is_valid_pdb_download(dest)) {
    if (file.exists(dest)) {
      unlink(dest)
    }
    return(FALSE)
  }

  TRUE
}

#' @keywords internal
is_valid_pdb_download = function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  size = file.info(path)$size
  if (is.na(size) || size <= 0L) {
    return(FALSE)
  }

  first_line = readLines(path, n = 1L, warn = FALSE)
  length(first_line) == 1L &&
    grepl(
      paste0(
        "^(HEADER|TITLE|COMPND|SOURCE|KEYWDS|EXPDTA|AUTHOR|",
        "REVDAT|JRNL|REMARK|MODEL |ATOM  |HETATM)"
      ),
      first_line
    )
}

#' @keywords internal
rcsb_pdb_destination = function(pdb_id, out_dir, filename = NULL) {
  if (is.null(filename)) {
    filename = paste0(tolower(pdb_id), ".pdb")
  } else {
    filename = validate_rcsb_filename(filename)
  }
  file.path(out_dir, filename)
}

#' @keywords internal
validate_rcsb_protein_query = function(protein) {
  if (!is.character(protein) || length(protein) != 1L || is.na(protein)) {
    stop("protein must be a single character value", call. = FALSE)
  }
  protein = trimws(protein)
  if (!nzchar(protein)) {
    stop("protein must not be empty", call. = FALSE)
  }
  protein
}

#' @keywords internal
validate_rcsb_filename = function(filename) {
  if (!is.character(filename) || length(filename) != 1L || is.na(filename)) {
    stop("filename must be NULL or a single character value", call. = FALSE)
  }
  if (!nzchar(filename)) {
    stop("filename must not be empty", call. = FALSE)
  }
  filename
}

#' @keywords internal
validate_rcsb_max_results = function(max_results) {
  if (
    length(max_results) != 1L ||
      is.na(max_results) ||
      !is.finite(max_results) ||
      max_results < 1L
  ) {
    stop("max_results must be a positive integer", call. = FALSE)
  }
  as.integer(max_results)
}

#' @keywords internal
validate_rcsb_flag = function(value, argument_name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be TRUE or FALSE", argument_name), call. = FALSE)
  }
  value
}

#' @keywords internal
is_rcsb_pdb_id = function(value) {
  grepl("^[0-9][A-Za-z0-9]{3}$", value)
}

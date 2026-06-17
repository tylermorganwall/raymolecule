# install.packages("jsonlite")

download_pdb_sample = function(n = 250, out_dir = "pdb_sample", seed = 1) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ids_url = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
  ids = jsonlite::fromJSON(ids_url)
  ids = toupper(ids)

  set.seed(seed)
  ids = sample(ids)

  files = character()
  failed = character()

  for (id in ids) {
    if (length(files) >= n) {
      break
    }

    dest = file.path(out_dir, paste0(tolower(id), ".pdb"))
    url = sprintf("https://files.rcsb.org/download/%s.pdb", id)

    ok = tryCatch(
      {
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)

        first_line = readLines(dest, n = 1, warn = FALSE)
        is_ok = file.exists(dest) &&
          file.info(dest)$size > 0 &&
          length(first_line) == 1 &&
          grepl(
            "^(HEADER|TITLE|COMPND|SOURCE|KEYWDS|EXPDTA|AUTHOR|REVDAT|JRNL|REMARK|MODEL |ATOM  |HETATM)",
            first_line
          )

        if (!is_ok) {
          unlink(dest)
        }

        is_ok
      },
      error = function(e) {
        if (file.exists(dest)) {
          unlink(dest)
        }
        FALSE
      }
    )

    if (ok) {
      files = c(files, dest)
    } else {
      failed = c(failed, id)
    }
  }

  list(files = files, failed = failed)
}

smoke_test_raymolecule = function(files, render = FALSE) {
  rows = lapply(files, function(path) {
    result = data.frame(
      file = path,
      parse_ok = FALSE,
      render_ok = FALSE,
      atoms = NA_integer_,
      residues = NA_integer_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )

    tryCatch(
      {
        model = raymolecule::read_pdb(path)

        result$parse_ok = TRUE
        result$atoms = if (!is.null(model$atoms)) nrow(model$atoms) else
          NA_integer_
        result$residues = if (!is.null(model$residues))
          nrow(model$residues) else NA_integer_

        has_protein_ca = !is.null(model$residues) &&
          "has_ca" %in% names(model$residues) &&
          any(model$residues$has_ca, na.rm = TRUE)

        if (has_protein_ca) {
          scene = raymolecule::generate_ribbon_scene(
            model,
            pathtrace = FALSE,
            show_waters = FALSE,
            cross_section_resolution = 8,
            subdivisions = 1
          )
        } else if (!is.null(model$atoms) && nrow(model$atoms) > 0L) {
          scene = raymolecule::generate_full_scene(
            model,
            pathtrace = FALSE
          )
        } else {
          stop("No renderable atoms were found in the PDB model")
        }

        result$render_ok = TRUE

        if (render) {
          raymolecule::render_model(scene, width = 256, height = 256)
        }

        result
      },
      error = function(e) {
        result$error = conditionMessage(e)
        result
      }
    )
  })

  do.call(rbind, rows)
}

# batch = download_pdb_sample(n = 10, out_dir = "pdb_sample", seed = 20260615)
# results = smoke_test_raymolecule(batch$files)

# table(results$parse_ok, results$render_ok, useNA = "ifany")
# head(results[!results$parse_ok | !results$render_ok, ], 20)

test_that("build_rcsb_protein_search_query builds an RCSB protein search", {
  query = build_rcsb_protein_search_query("hemoglobin", max_results = 5L)

  expect_equal(query$return_type, "entry")
  expect_equal(query$query$type, "group")
  expect_equal(query$query$logical_operator, "and")
  expect_equal(query$query$nodes[[1]]$service, "full_text")
  expect_equal(query$query$nodes[[1]]$parameters$value, "hemoglobin")
  expect_equal(query$query$nodes[[2]]$service, "text")
  expect_equal(
    query$query$nodes[[2]]$parameters$attribute,
    "rcsb_entry_info.polymer_entity_count_protein"
  )
  expect_equal(query$query$nodes[[2]]$parameters$operator, "greater")
  expect_equal(query$request_options$results_content_type, list("experimental"))
  expect_equal(query$request_options$paginate$rows, 5L)
})

test_that("extract_rcsb_pdb_ids handles common RCSB response shapes", {
  expect_equal(
    extract_rcsb_pdb_ids(list(
      result_set = list(
        list(identifier = "4hhb", score = 1),
        list(identifier = "1abc", score = 0.5)
      )
    )),
    c("4HHB", "1ABC")
  )

  expect_equal(
    extract_rcsb_pdb_ids(list(result_set = c("4hhb", "1abc"))),
    c("4HHB", "1ABC")
  )

  expect_equal(
    extract_rcsb_pdb_ids(list(
      result_set = data.frame(
        identifier = c("4hhb", NA, ""),
        stringsAsFactors = FALSE
      )
    )),
    "4HHB"
  )

  expect_equal(extract_rcsb_pdb_ids(list()), character())
})

test_that("is_valid_pdb_download validates downloaded PDB-like files", {
  file = tempfile(fileext = ".pdb")
  on.exit(unlink(file), add = TRUE)

  writeLines("HEADER    TEST STRUCTURE", file)
  expect_true(is_valid_pdb_download(file))

  writeLines("404: Not Found", file)
  expect_false(is_valid_pdb_download(file))

  unlink(file)
  expect_false(is_valid_pdb_download(file))
})

test_that("download_pdb downloads direct PDB IDs and annotates the path", {
  temp_dir = tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  testthat::local_mocked_bindings(
    download_rcsb_pdb_file = function(pdb_id, dest, overwrite) {
      expect_equal(pdb_id, "4HHB")
      expect_false(overwrite)
      writeLines("HEADER    MOCK PDB", dest)
      TRUE
    }
  )

  path = download_pdb("4hhb", out_dir = temp_dir)

  expect_equal(basename(path), "4hhb.pdb")
  expect_true(file.exists(path))
  expect_equal(attr(path, "pdb_id"), "4HHB")
  expect_equal(attr(path, "query"), "4hhb")
  expect_equal(attr(path, "matched_ids"), "4HHB")
})

test_that("download_pdb searches by protein name and tries later valid hits", {
  temp_dir = tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  attempts = character()

  testthat::local_mocked_bindings(
    rcsb_search_pdb_ids = function(protein, max_results) {
      expect_equal(protein, "hemoglobin")
      expect_equal(max_results, 2L)
      c("1BAD", "4HHB")
    },
    download_rcsb_pdb_file = function(pdb_id, dest, overwrite) {
      attempts <<- c(attempts, pdb_id)
      if (pdb_id == "1BAD") {
        return(FALSE)
      }
      writeLines("HEADER    MOCK PDB", dest)
      TRUE
    }
  )

  path = download_pdb("hemoglobin", out_dir = temp_dir, max_results = 2L)

  expect_equal(attempts, c("1BAD", "4HHB"))
  expect_equal(basename(path), "4hhb.pdb")
  expect_equal(attr(path, "pdb_id"), "4HHB")
  expect_equal(attr(path, "matched_ids"), c("1BAD", "4HHB"))
})

test_that("download_pdb validates inputs before searching", {
  expect_error(download_pdb(character()), "protein must be a single")
  expect_error(download_pdb(""), "protein must not be empty")
  expect_error(download_pdb("hemoglobin", max_results = 0), "max_results")
  expect_error(download_pdb("hemoglobin", verbose = NA), "verbose")
  expect_error(download_pdb("hemoglobin", overwrite = NA), "overwrite")
})

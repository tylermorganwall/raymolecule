#' Get Molecule
#'
#' Loads the structure of a molecule by fetching an SDF file from Pubchem, which can be piped to generate_full_scene
#'
#' @param molecule A character variable of a compound name or a numeric variable of an official compound ID
#'
#' @return List giving the atom locations and the connections between atoms.
#' @export
#'
#' @examples
#' get_molecule("caffeine")
#' get_molecule(5757) #estradiol (aka estrogen)
#' get_molecule("aspirin")
get_molecule = function(molecule) {
  #molecule details
  if (is.numeric(molecule)) {    #assume it is a CID
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", molecule, "/SDF")
  } else if (is.character(molecule)) {    #assume it is a compound name
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", molecule, "/SDF")
  } else {
    stop("Provide a quoted compound name (character) or an unquoted compound ID (numeric)")
  }

  status <- httr::GET(url) %>%
    httr::status_code()

  if (status == 200) {
    molecule_sdf <-
      raymolecule::read_sdf(url)
  } else {
    stop("Cannot find your compound. Check the compound name or compound ID (CID)")
  }
  return(molecule_sdf)
}

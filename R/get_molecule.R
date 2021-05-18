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
#' \donttest{
#' get_molecule("caffeine")
#' get_molecule(5757) #estradiol (aka estrogen)
#' get_molecule("testosterone")
#' get_molecule("aspirin")
#' }
get_molecule = function(molecule) {
  #molecule details
  if (is.numeric(molecule)) {    #assume it is a CID
    url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF",molecule)
  } else if (is.character(molecule)) {    #assume it is a compound name
    url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/SDF",molecule)
  } else {
    stop("Provide a quoted compound name (character) or an unquoted compound ID (numeric)")
  }

  status = httr::status_code(httr::GET(url))

  if (status == 200) {
    molecule_sdf = read_sdf(url)
  } else {
    stop(sprintf("Cannot find compound '%s'. Check the compound name or compound ID (CID)", molecule))
  }
  return(molecule_sdf)
}

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
#' get_molecule("caffeine") %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' #estradiol (aka estrogen)
#' get_molecule(5757) %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' get_molecule("testosterone") %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' get_molecule("aspirin") %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' get_molecule("rutoside") %>%
#'   generate_full_scene() %>%
#'   render_model()
#'
#' #If the 3D SDF doesn't exist, this function will pull the 2D SDF and inform the user
#' get_molecule("cyanocobalamin") %>%
#'   generate_full_scene() %>%
#'   render_model()
#' }
get_molecule = function(molecule) {
  #Check for 3D model first
  if (is.numeric(molecule)) {    #assume it is a CID
    url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF/?record_type=3d",molecule)
  } else if (is.character(molecule)) {    #assume it is a compound name
    url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/SDF/?record_type=3d",molecule)
  } else {
    stop("Provide a quoted compound name (character) or an unquoted compound ID (numeric)")
  }

  status = httr::status_code(httr::GET(url))

  #Check for 2D model if 3D doesn't return SDF
  if(status != 200) {
    if (is.numeric(molecule)) {    #assume it is a CID
      url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF",molecule)
    } else if (is.character(molecule)) {    #assume it is a compound name
      url = sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/SDF",molecule)
    } else {
      stop("Provide a quoted compound name (character) or an unquoted compound ID (numeric)")
    }
  }
  status = httr::status_code(httr::GET(url))

  if (status == 200) {
    molecule_sdf = read_sdf(url)
  } else {
    stop(sprintf("Cannot find compound '%s'. Check the compound name or compound ID (CID)", molecule))
  }
  return(molecule_sdf)
}

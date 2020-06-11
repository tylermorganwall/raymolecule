#' Get Example Molecule
#'
#' Loads the structure of the built-in molecules. All SDF files obtained from Pubchem.
#' txt extension only included to pass R CHECK.
#'
#' @param molecule One of the built-in SDF files. These are "benzene", "buckyball",
#' "caffeine", "capsaicin", "cinnemaldehyde", "geraniol", "luciferin", "morphine",
#' "penicillin", "pfoa", "skatole", "tubocurarine_chloride".
#'
#' @return List giving the atom locations and the connections between atoms.
#' @export
#'
#' @examples
#' get_example_molecule("benzene")
#' get_example_molecule("cinnemaldehyde")
#' get_example_molecule("geraniol")
get_example_molecule = function(molecule) {
  if(!molecule %in% c("benzene", "buckyball","caffeine", "capsaicin", "cinnemaldehyde", "geraniol", "luciferin", "morphine",
                     "penicillin", "pfoa", "skatole", "tubocurarine_chloride")) {
    stop("molecule ", molecule, "not found")
  }
  system.file("extdata", sprintf("%s.txt",molecule), package="raymolecule")
}

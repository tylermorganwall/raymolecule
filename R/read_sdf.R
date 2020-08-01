#' Read SDF File
#'
#' Reads an SDF file and extracts the 3D molecule model
#'
#' @param filename Filename to the sdf file.
#'
#' @return List giving the atom locations and the connections between atoms.
#' @export
#'
#' @examples
#' #This assumes a hypothetical SDF file in your working directory:
#' \dontrun{
#' read_pdb("molecule.sdf") %>%
#'   generate_full_scene() %>%
#'   render_model()
#' }
read_sdf = function(filename) {
  as.numeric = function(...) {
    suppressWarnings(base::as.numeric(...))
  }
  con = file(filename,"rt")
  on.exit(close(con))
  readLines(con,n=3)
  values = as.numeric(unlist(strsplit(trimws(readLines(con,n=1)), split = "\\s", perl=TRUE)))
  molecules = values[1]
  bonds = values[2]
  moleculelist = strsplit(trimws(readLines(con,n=molecules)), split = "\\s+", perl=TRUE)
  moleculedf = as.data.frame(do.call(rbind,moleculelist), stringsAsFactors = FALSE)
  moleculedf = moleculedf[,1:4]
  moleculedf= cbind(moleculedf,list(index=1:nrow(moleculedf)))
  moleculedf[,1] = as.numeric(moleculedf[,1])
  moleculedf[,2] = as.numeric(moleculedf[,2])
  moleculedf[,3] = as.numeric(moleculedf[,3])
  colnames(moleculedf) = c("x","y","z","type","index")
  bondlist = strsplit(trimws(readLines(con,n=bonds)), split = "\\s+", perl=TRUE)
  bonddf = as.data.frame(do.call(rbind,bondlist), stringsAsFactors = FALSE)
  bonddf[,1] = as.numeric(bonddf[,1])
  bonddf[,2] = as.numeric(bonddf[,2])
  bonddf[,3] = as.numeric(bonddf[,3])
  colnames(bonddf) = c("from", "to", "number")
  return(list(atoms = moleculedf, bonds = bonddf))
}

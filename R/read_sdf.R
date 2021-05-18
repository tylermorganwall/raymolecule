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
  con = file(filename,"rb")
  on.exit(close(con))
  readLines(con,n=3)
  molecules = as.numeric(readChar(con,nchars = 3))
  bonds = as.numeric(readChar(con,nchars = 3))
  readLines(con,n=1)

  moleculelist = list()
  for(i in seq_len(molecules)) {
    moleculelist[[i]] = matrix(trimws(c(readChar(con,nchars = 10),
                                 readChar(con,nchars = 10),
                                 readChar(con,nchars = 10),
                                 readChar(con,nchars = 3),
                                 readChar(con,nchars = 3))),ncol=5)
    readLines(con,n=1)
  }

  moleculedf = as.data.frame(do.call(rbind,moleculelist), stringsAsFactors = FALSE)
  moleculedf = moleculedf[,1:4]
  moleculedf= cbind(moleculedf,list(index=1:nrow(moleculedf)))
  moleculedf[,1] = as.numeric(moleculedf[,1])
  moleculedf[,2] = as.numeric(moleculedf[,2])
  moleculedf[,3] = as.numeric(moleculedf[,3])
  colnames(moleculedf) = c("x","y","z","type","index")
  bondlist = list()
  for(i in seq_len(bonds)) {
    bondlist[[i]] = matrix(as.numeric(trimws(c(readChar(con,nchars = 3),
                                        readChar(con,nchars = 3),
                                        readChar(con,nchars = 3)))),ncol=3)
    readLines(con,n=1)
  }
  bonddf = as.data.frame(do.call(rbind,bondlist), stringsAsFactors = FALSE)
  colnames(bonddf) = c("from", "to", "number")
  return(list(atoms = moleculedf, bonds = bonddf))
}

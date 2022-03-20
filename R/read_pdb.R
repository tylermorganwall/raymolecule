#' Read PDB File
#'
#' Reads an PDB file and extracts the atom locations and bonds (does not include any other structual
#' information currently). This pulls out ATOM and HETAHM records by default, along with available
#' connections.
#'
#' @param filename Path to the PDB file.
#' @param atom Default `TRUE`. Whether to pull out standard residue (ATOM) records.
#' @param nsr Default `TRUE`. Whether to pull out nonstandard residue (HETAHM) records.
#'
#' @return List giving the atom locations.
#' @export
#'
#' @examples
#' #This assumes a hypothetical PDB file in your working directory:
#' if(file.exists("3nir.pdb")) {
#'   read_pdb("3nir.pdb") %>%
#'     generate_full_scene() %>%
#'     render_model()
#' }
read_pdb = function(filename, atom = TRUE, nsr = TRUE) {
  as.numeric = function(...) {
    suppressWarnings(base::as.numeric(...))
  }
  con = file(filename,"rb")
  on.exit(close(con))
  type = ""
  while(!trimws(type[[1]]) %in% c("ATOM", "HETATM")) {
    type = trimws(readChar(con,nchars=6))
    if(!type %in% c("ATOM", "HETATM")) {
      junk = readLines(con,n=1)
    }
  }
  number = as.integer(readChar(con,nchars=5))
  other = readChar(con,nchars=19)
  x = as.numeric(readChar(con,nchars=8))
  y = as.numeric(readChar(con,nchars=8))
  z = as.numeric(readChar(con,nchars=8))
  other = readChar(con,nchars=12)
  type = trimws(readChar(con,nchars=12))
  junk = readLines(con,n=1)
  xlist = list()
  ylist = list()
  zlist = list()
  typelist = list()
  indexlist = list()
  bondlist = list()
  xlist[[1]] = x
  ylist[[1]] = y
  zlist[[1]] = z
  typelist[[1]] = type
  indexlist[[1]] = number
  counter = 2
  counter_bond = 1
  val = ""
  while("END" != val) {
    val = trimws(readChar(con,nchars=6))
    if(("ATOM" == val && atom) || ("HETATM" == val && nsr)) {
      number = as.integer(readChar(con,nchars=5))
      other = readChar(con,nchars=19)
      x = as.numeric(readChar(con,nchars=8))
      y = as.numeric(readChar(con,nchars=8))
      z = as.numeric(readChar(con,nchars=8))
      other = readChar(con,nchars=12+10)
      type = trimws(readChar(con,nchars=2))
      xlist[[counter]] = x
      ylist[[counter]] = y
      zlist[[counter]] = z
      typelist[[counter]] = type
      indexlist[[counter]] = number
      counter = counter + 1
    }
    if("CONECT" == val) {
      source = as.numeric(readChar(con,nchars=5))
      dest = as.numeric(c(readChar(con,nchars=5),readChar(con,nchars=5),readChar(con,nchars=5)))
      dest = dest[!is.na(dest)]
      if(length(dest) > 0) {
        bondlist[[counter_bond]] = data.frame(from = source, to = dest, number = 1)
        counter_bond = counter_bond + 1
      }
    }
    junk = readLines(con,n=1)
  }
  moleculedf = data.frame(x=unlist(xlist),y=unlist(ylist), z=unlist(zlist),
                          type=unlist(typelist),index = unlist(indexlist))
  moleculedf = moleculedf[!is.na(moleculedf$x),]
  bonddf = do.call(rbind,bondlist)
  return(list(atoms=moleculedf, bonds = bonddf))
}

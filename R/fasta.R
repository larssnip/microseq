#' @name readFasta
#' @title Read and write FASTA files
#' @aliases readFasta writeFasta
#' 
#' @description Reads and writes biological sequences (DNA, RNA, protein) in the FASTA format.
#' 
#' @usage readFasta(in.file)
#' writeFasta(fdta, out.file, width = 80)
#' 
#' @param in.file url/directory/name of FASTA file to read.
#' @param fdta A \code{\link{tibble}} with sequence data, see \sQuote{Details} below.
#' @param out.file Name of FASTA file to create.
#' @param width Number of characters per line, or 0 for no linebreaks.
#' 
#' @details These functions handle input/output of sequences in the commonly used FASTA format.
#' For every sequence it is presumed there is one Header-line starting with a \sQuote{>}.
#' 
#' The sequences are stored in a \code{\link{tibble}}, opening up all the possibilities in R for
#' fast and easy manipulations. The content of the file is stored as two columns, \samp{Header}
#' and \samp{Sequence}. If other columns are added, these will be ignored by \code{\link{writeFasta}}.
#' 
#' Setting \code{width = 0} in \code{\link{writeFasta}} results in no linebreaks in the sequences
#' (one sequence per line).
#' 
#' @return \code{\link{readFasta}} returns a \code{\link{tibble}} with the contents of the FASTA
#' file stored in two columns of text. The first, named \samp{Header}, contains
#' the headerlines and the second, named \samp{Sequence}, contains the sequences.
#' 
#' \code{\link{writeFasta}} produces a FASTA file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readFastq}}.
#' 
#' @examples
#' \dontrun{
#' # We need a FASTA-file to read, here is one example file:
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.ffn")
#' 
#' # Read and write
#' fdta <- readFasta(fa.file)
#' ok <- writeFasta(fdta[4:5,], out.file = "delete_me.fasta")
#' 
#' # Make use of dplyr to copy parts of the file to another file
#' readFasta(fa.file) %>% 
#'   filter(str_detect(Sequence, "TGA$")) %>% 
#'   writeFasta(out.file = "TGAstop.fasta", width = 0) -> ok
#' }
#' 
#' @keywords sequence FASTA
#' 
#' @importFrom tibble as_tibble
#' @importFrom stringr str_c
#' @importFrom dplyr %>% 
#' 
#' @export readFasta
#' @export writeFasta
#' 
readFasta <- function(in.file){
  if(!is.character(in.file) | length(in.file)>1) stop("The argument in.file must be a single text (a filename)")
  if(file.exists(in.file)){
    in.file <- normalizePath(in.file)
    read_fasta(in.file) %>% 
      as_tibble() -> fdta
    return(fdta)
  } else {
    stop("Cannot find ", in.file, ", please correct path and/or file name")
  }
}
writeFasta <- function(fdta, out.file, width = 80){
  if(!("Header" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Header is lacking\n")
  }
  if(!("Sequence" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Sequence is lacking\n")
  }
  out.file <- file.path(normalizePath(dirname(out.file)),
                        basename(out.file))
  if(width == 0){
    str_c(">", fdta$Header) %>% 
      cbind(fdta$Sequence) %>% 
      t() %>% 
      as.vector() %>% 
      writeLines(con = out.file)
    status <- TRUE
  } else {
    status <- write_fasta(fdta$Header, fdta$Sequence, out.file, width)
  }
  invisible(status)
}

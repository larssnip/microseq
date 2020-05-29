#' @name readFasta
#' @title Read and write FASTA files
#' @aliases readFasta writeFasta
#' 
#' @description Reads and writes biological sequences (DNA, RNA, protein) in the FASTA format.
#' 
#' @usage readFasta(in.file)
#' writeFasta(fdta, out.file, width = 0)
#' 
#' @param in.file url/directory/name of (gzipped) FASTA file to read.
#' @param fdta A \code{\link{tibble}} with sequence data, see \sQuote{Details} below.
#' @param out.file Name of (gzipped) FASTA file to create.
#' @param width Number of characters per line, or 0 for no linebreaks.
#' 
#' @details These functions handle input/output of sequences in the commonly used FASTA format.
#' For every sequence it is presumed there is one Header-line starting with a \sQuote{>}. If
#' filenames (\code{in.file} or \code{out.file}) have the extension \code{.gz} they will automatically be
#' compressed/uncompressed.
#' 
#' The sequences are stored in a \code{\link{tibble}}, opening up all the possibilities in R for
#' fast and easy manipulations. The content of the file is stored as two columns, \samp{Header}
#' and \samp{Sequence}. If other columns are added, these will be ignored by \code{\link{writeFasta}}.
#' 
#' The default \code{width = 0} in \code{\link{writeFasta}} results in no linebreaks in the sequences
#' (one sequence per line).
#' 
#' @return \code{\link{readFasta}} returns a \code{\link{tibble}} with the contents of the (gzipped) FASTA
#' file stored in two columns of text. The first, named \samp{Header}, contains
#' the headerlines and the second, named \samp{Sequence}, contains the sequences.
#' 
#' \code{\link{writeFasta}} produces a (gzipped) FASTA file.
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
#'   writeFasta(out.file = "TGAstop.fasta", width = 80) -> ok
#' }
#' 
#' @keywords sequence FASTA
#' 
#' @importFrom tibble as_tibble
#' @importFrom stringr str_c
#' @importFrom dplyr %>% 
#' @importFrom data.table fread fwrite
#' 
#' @export readFasta
#' @export writeFasta
#' 
readFasta <- function(in.file){
  if(!is.character(in.file) | length(in.file) > 1) stop("The argument in.file must be a single text (a filename)")
  if(file.exists(in.file)){
    in.file <- normalizePath(in.file)
    tbl <- fread(in.file, header = F, sep = "\t", data.table = F)
    hh <- grep("^>", tbl[,1])
    hhh <- c(hh, nrow(tbl) + 1)
    tibble(Header = tbl[hh,1],
           Sequence = sapply(1:length(hh), function(i){str_c(tbl[(hhh[i]+1):(hhh[i+1]-1),1], collapse = "")})) -> fdta
    return(fdta)
  } else {
    stop("Cannot find ", in.file, ", please correct path and/or file name")
  }
}
writeFasta <- function(fdta, out.file, width = 0){
  if(!("Header" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Header is lacking\n")
  }
  if(!("Sequence" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Sequence is lacking\n")
  }
  out.file <- file.path(normalizePath(dirname(out.file)),
                        basename(out.file))
  lst <- vector("list", 2 * nrow(fdta))
  lst[seq(1, 2 * nrow(fdta), 2)] <- as.list(str_c(">", fdta$Header))
  if(width > 0){
    lst[seq(2, 2 * nrow(fdta), 2)] <- lapply(fdta$Sequence, chop, width)
  } else {
    lst[seq(2, 2 * nrow(fdta), 2)] <- as.list(fdta$Sequence)
  }
  fwrite(tibble(unlist(lst)), file = out.file, compress = "auto", col.names = F)
}

chop <- function(sequence, width){
  n <- str_length(sequence)
  s1 <- seq(1, n, width)
  s2 <- pmin(s1 + width - 1, n)
  return(str_sub(sequence, s1, s2))
}

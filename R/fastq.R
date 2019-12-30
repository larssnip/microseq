#' @name readFastq
#' @title Read and write FASTQ files
#' @aliases readFastq writeFastq
#' 
#' @description Reads and writes files in the FASTQ format.
#' 
#' @usage readFastq(in.file)
#' writeFastq(fdta, out.file)
#' 
#' @param in.file url/directory/name of FASTQ file to read, possibly compressed to .gz.
#' @param fdta FASTQ object to write.
#' @param out.file url/directory/name of FASTQ file to write.
#' 
#' @details These functions handle input/output of sequences in the commonly used FASTQ format,
#' typically used for storing DNA sequences (reads) after sequencing. 
#' 
#' The sequences are stored in a \code{\link{tibble}}, opening up all the possibilities in R for
#' fast and easy manipulations. The content of the file is stored as three columns, \samp{Header},
#' \samp{Sequence} and \samp{Quality}. If other columns are added, these will be ignored by
#' \code{\link{writeFastq}}.
#' 
#' @note These functions will only handle files where each entry spans one single line, i.e. not the
#' (uncommon) multiline FASTQ format.
#' 
#' @return \code{\link{readFastq}} returns a \code{\link{tibble}} with the contents of the FASTQ
#' file stored in three columns of text. The first, named \samp{Header}, contains
#' the headerlines, the second, named \samp{Sequence}, contains the sequences and the third, named 
#' \samp{Quality} contains the base quality scores.
#' 
#' \code{\link{writeFastq}} produces an uncompressed FASTQ file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso code{\link{readFasta}}.
#' 
#' @examples 
#' \dontrun{
#' # We need a FASTQ-file to read, here is one example file:
#' fq.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fastq")
#' 
#' # Read and write
#' fdta <- readFastq(fq.file)
#' ok <- writeFastq(fdta[1:3,], out.file = "delete_me.fq")
#' 
#' # Make use of dplyr to copy parts of the file to another file
#' readFastq(fq.file) %>% 
#'   mutate(Length = str_length(Sequence)) %>% 
#'   filter(Length > 200) %>% 
#'   writeFasta(out.file = "long_reads.fa") # writing to FASTA file
#' }
#' 
#' @keywords sequence FASTQ
#' 
#' @useDynLib microseq, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr %>%
#' @importFrom stringr str_detect
#' 
#' @export readFastq
#' @export writeFastq
#' 
readFastq <- function(in.file){
  if(!is.character(in.file) | length(in.file)>1) stop("The argument in.file must be a single text (a filename)")
  if(file.exists(in.file)){
    in.file <- normalizePath(in.file)
    if(str_detect(in.file, "gz$")){
      con <- gzfile(in.file, open = "rt")
      lines <- readLines(con)
      close(con)
      fq <- tibble(Header = lines[seq(1, length(lines), 4)],
                   Sequence = lines[seq(2, length(lines), 4)],
                   Quality = lines[seq(4, length(lines), 4)])
    } else {
      fq <- read_fastq(in.file)
    }
    return(as_tibble(fq))
  } else {
    stop("Cannot find ", in.file, ", please correct path and/or file name")
  }
}
writeFastq <- function(fdta, out.file){
  if(!("Header" %in% colnames(fdta))){
    stop("This is not a fastq object, column Header is lacking\n")
  }
  if(!("Sequence" %in% colnames(fdta))){
    stop("This is not a fastq object, column Sequence is lacking\n")
  }
  if(!("Quality" %in% colnames(fdta))){
    stop("This is not a fastq object, column Quality is lacking\n")
  }
  out.file <- file.path(normalizePath(dirname(out.file)),
                        basename(out.file))
  status <- write_fastq(fdta$Header, fdta$Sequence, fdta$Quality, out.file)
  invisible(status)
}




















#' #' @name plot.Fastq
#' #' @title Plotting and summary of \code{Fastq} objects
#' #' @aliases plot.Fastq summary.Fastq
#' #' 
#' #' @description Generic functions for plotting and printing the content of a \code{Fastq} object.
#' #' 
#' #' @param x A \code{Fastq} object, see below.
#' #' @param y not used.
#' #' @param object A \code{Fastq} object, see below.
#' #' @param col Color of bar interiors.
#' #' @param border Color of bar borders.
#' #' @param ... Optional graphical arguments.
#' #' 
#' #' @details  A \code{Fastq} object contains biological sequences in the FASTQ format. It is a small (S3)
#' #' extension to a \code{data.frame}. It is actually a \code{data.frame} containing at least three text columns
#' #' named \samp{Header}, \samp{Sequence} and \samp{Quality}. The \samp{Header} column contains the headerlines
#' #' for each sequence, the \samp{Sequence} columns the sequences themselves and the \samp{Quality} the base-quality
#' #' sequence. A \code{Fastq} object is typically created by reading
#' #' a FASTQ formatted file into R by \code{\link{readFastq}}.
#' #' 
#' #' A \code{Fastq} object can be treated as a \code{data.frame}, which makes it quick and easy to search both
#' #' for specific regular expressions, sort or re-arrange the ordering of the sequences,
#' #' extract subsets or add new data to an existing \code{Fastq} object.
#' #' 
#' #' A \code{Fastq} object can also be treated as a \code{Fasta} object. Using \code{\link{writeFasta}}
#' #' will write the \code{Fastq} to a file in \code{Fasta} format.
#' #' 
#' #' The \code{plot.Fastq} function will display the content of the \code{Fastq} object as a histogram over
#' #' the lengths of the sequences. 
#' #' 
#' #' The \code{summary.Fastq} function will display a text giving the number of sequences and the alphabet,
#' #' i.e. listing all unique symbols found in the file.
#' #' 
#' #' @author Lars Snipen and Kristian Hovde Liland.
#' #' 
#' #' @seealso \code{\link{readFastq}}, \code{\link{writeFastq}}.
#' #' 
#' #' @examples # See the examples in the Help-file for readFastq/writeFastq
#' #' 
#' #' @importFrom graphics plot hist
#' #' 
#' #' @method plot Fastq
#' #' @export
#' plot.Fastq <- function(x, y = NULL, col = "tan4", border = "tan4", ...){
#'   nc <- nchar(x$Sequence)
#'   hist(nc, breaks = length(nc)^(1/3), col = col, xlab = "Sequence length",
#'         ylab = "Number of sequences", main = "Fastq object" )
#' }
#' 
#' #' @rdname plot.Fastq
#' #' @method summary Fastq
#' #' @export
#' summary.Fastq <- function(object, ...){
#'   n.seq <- nrow(object)
#'   alphabet <- unique(unlist(strsplit(object$Sequence, split = "")))
#'   cat("Fastq formatted sequence data containing", n.seq, "sequences\n")
#'   cat("Alphabet:", sort(alphabet), "\n")
#' }
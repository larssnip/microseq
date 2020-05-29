#' @name readFastq
#' @title Read and write FASTQ files
#' @aliases readFastq writeFastq
#' 
#' @description Reads and writes files in the FASTQ format.
#' 
#' @usage readFastq(in.file)
#' writeFastq(fdta, out.file)
#' 
#' @param in.file url/directory/name of (gzipped) FASTQ file to read.
#' @param fdta FASTQ object to write.
#' @param out.file url/directory/name of (gzipped) FASTQ file to write.
#' 
#' @details These functions handle input/output of sequences in the commonly used FASTQ format,
#' typically used for storing DNA sequences (reads) after sequencing. If
#' filenames (\code{in.file} or \code{out.file}) have the extension \code{.gz} they will automatically be
#' compressed/uncompressed. 
#' 
#' The sequences are stored in a \code{\link{tibble}}, opening up all the possibilities in R for
#' fast and easy manipulations. The content of the file is stored as three columns, \samp{Header},
#' \samp{Sequence} and \samp{Quality}. If other columns are added, these will be ignored by
#' \code{\link{writeFastq}}.
#' 
#' @note These functions will only handle files where each entry spans one single line, i.e. not the
#' (uncommon) multiline FASTQ format.
#' 
#' @return \code{\link{readFastq}} returns a \code{\link{tibble}} with the contents of the (gzipped) FASTQ
#' file stored in three columns of text. The first, named \samp{Header}, contains
#' the headerlines, the second, named \samp{Sequence}, contains the sequences and the third, named 
#' \samp{Quality} contains the base quality scores.
#' 
#' \code{\link{writeFastq}} produces a (gzipped) FASTQ file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso code{\link{readFasta}}.
#' 
#' @examples 
#' \dontrun{
#' # We need a FASTQ-file to read, here is one example file:
#' fq.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fastq.gz")
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
#' @importFrom data.table fread fwrite
#' @importFrom rlang .data
#' 
#' @export readFastq
#' @export writeFastq
#' 
readFastq <- function(in.file){
  if(!is.character(in.file) | length(in.file)>1) stop("The argument in.file must be a single text (a filename)")
  if(file.exists(in.file)){
    in.file <- normalizePath(in.file)
    tbl <- fread(in.file, header = F, sep = "\t", data.table = F)
    tibble(Header   = str_remove(tbl[seq(1, nrow(tbl), 4),1], "^@"),
           Sequence = tbl[seq(2, nrow(tbl), 4),1],
           Quality  = tbl[seq(4, nrow(tbl), 4),1]) -> fdta
    return(fdta)
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
  fdta %>% 
    mutate(Plus = "+") %>% 
    mutate(Header = str_c("@", .data$Header)) %>% 
    select(.data$Header, .data$Sequence, .data$Plus, .data$Quality) %>% 
    as.matrix() %>% 
    t() %>% 
    as.character() %>% 
    as_tibble() -> tbl
  fwrite(tbl, file = out.file, compress = "auto", col.names = F)
}

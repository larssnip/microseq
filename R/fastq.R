#' @name readFastq
#' @title Read and write FASTQ files
#' @aliases readFastq writeFastq
#' 
#' @description Reads and writes files in the FASTQ format.
#' 
#' @usage readFastq(in.file, Sanger = FALSE)
#' writeFastq(fdta, out.file)
#' 
#' @param in.file url/directory/name of FASTQ file to read.
#' @param Sanger logical indicating if old, multi-line Sanger format is used (default = FALSE).
#' @param fdta FASTQ object to write.
#' @param out.file url/directory/name of FASTQ file to write.
#' 
#' @details These functions handle input/output of sequences in the commonly used FASTQ format,
#' typically used for storing DNA sequences (reads) after sequencing. 
#' 
#' The sequences are stored in a \code{Fastq} object. This is an extension of a \code{data.frame}
#' containing three text-columns named \samp{Header}, \samp{Sequence} and \samp{Quality}. If
#' other columns are present, these will be ignored by \code{writeFastq}.
#' 
#' The \code{Fastq} object can be treated as a \code{data.frame}, but the generic functions
#' \code{\link{plot.Fastq}} and \code{\link{summary.Fastq}} are defined. The \code{data.frame}
#' property makes it straightforward to manipulate all headers or all sequences, or to extract
#' or delete entries (rows), or to merge several data sets using \code{\link{rbind}}.
#' 
#' A \code{Fastq} object can also be treated as a \code{Fasta} object. Using \code{\link{writeFasta}}
#' will write the \code{Fastq} object to a file in \code{Fasta} format.
#' 
#' @return \code{\link{readFastq}} returns a \code{Fastq} object with the contents of the FASTQ file.
#' This is an extension to a \code{data.frame} and contains three columns of text. The first,
#' named \samp{Header}, contains the headerlines and the second, named \samp{Sequence}, contains
#' the sequences and the third, named \samp{Quality}, contains the base quality scores.
#' 
#' @return \code{\link{readFastq}} returns a \code{Fastq} object with the header, sequence and
#' quality part of the FASTQ file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso code{\link{readFasta}}, \code{\link{plot.Fastq}}, \code{\link{summary.Fastq}}.
#' 
#' @examples 
#' \dontrun{
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fastq")
#' fdta <- readFastq(ex.file)
#' summary(fdta)
#' }
#' 
#' @keywords sequence FASTQ Fastq
#' 
#' @useDynLib microseq
#' @importFrom Rcpp evalCpp
#' 
#' @export readFastq
#' @export writeFastq
#' 
readFastq <- function( in.file, Sanger = FALSE ){
  if( file.exists( in.file ) ){
    in.file <- normalizePath( in.file )
    if(Sanger){
      fq <- read_fastq_Sanger( in.file )
    } else {
      fq <- read_fastq( in.file )
    }
    fq <- as.data.frame( fq, stringsAsFactors = FALSE )
    class( fq ) <- c( "Fastq", "data.frame" )
    return( fq )
  } else {
    stop( "Cannot find ", in.file, ", please correct path and/or file name" )
  }
}

writeFastq <- function( fdta, out.file ){
  cn <- colnames( fdta )
  if( !("Header" %in% cn) || !("Sequence" %in% cn) || !("Quality" %in% cn) ){
    stop( "This is not a Fastq object, Header, Sequence or Quality is lacking\n" )
  }
  status <- write_fastq( fdta$Header, fdta$Sequence, fdta$Quality, out.file )
  invisible( status )
}


#' @name plot.Fastq
#' @title Plotting and summary of \code{Fastq} objects
#' @aliases plot.Fastq summary.Fastq
#' 
#' @description Generic functions for plotting and printing the content of a \code{Fastq} object.
#' 
#' @param x A \code{Fastq} object, see below.
#' @param y not used.
#' @param object A \code{Fastq} object, see below.
#' @param col Color of bar interiors.
#' @param border Color of bar borders.
#' @param ... Optional graphical arguments.
#' 
#' @details  A \code{Fastq} object contains biological sequences in the FASTQ format. It is a small (S3)
#' extension to a \code{data.frame}. It is actually a \code{data.frame} containing at least three text columns
#' named \samp{Header}, \samp{Sequence} and \samp{Quality}. The \samp{Header} column contains the headerlines
#' for each sequence, the \samp{Sequence} columns the sequences themselves and the \samp{Quality} the base-quality
#' sequence. A \code{Fastq} object is typically created by reading
#' a FASTQ formatted file into R by \code{\link{readFastq}}.
#' 
#' A \code{Fastq} object can be treated as a \code{data.frame}, which makes it quick and easy to search both
#' for specific regular expressions, sort or re-arrange the ordering of the sequences,
#' extract subsets or add new data to an existing \code{Fastq} object.
#' 
#' A \code{Fastq} object can also be treated as a \code{Fasta} object. Using \code{\link{writeFasta}}
#' will write the \code{Fastq} to a file in \code{Fasta} format.
#' 
#' The \code{plot.Fastq} function will display the content of the \code{Fastq} object as a histogram over
#' the lengths of the sequences. 
#' 
#' The \code{summary.Fastq} function will display a text giving the number of sequences and the alphabet,
#' i.e. listing all unique symbols found in the file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readFastq}}, \code{\link{writeFastq}}.
#' 
#' @examples # See the examples in the Help-file for readFastq/writeFastq
#' 
#' @importFrom graphics plot hist
#' 
#' @method plot Fastq
#' @export
plot.Fastq <- function( x, y=NULL, col="tan4", border="tan4", ... ){
  nc <- nchar( x$Sequence )
  hist( nc, breaks=length( nc )^(1/3), col="tan4", xlab="Sequence length",
        ylab="Number of sequences", main="Fastq object" )
}

#' @rdname plot.Fastq
#' @method summary Fastq
#' @export
summary.Fastq <- function( object, ... ){
  n.seq <- nrow( object )
  alphabet <- unique( unlist( strsplit( object$Sequence, split="" ) ) )
  cat( "Fastq formatted sequence data containing", n.seq, "sequences\n" )
  cat( "Alphabet:", sort( alphabet ), "\n" )
}
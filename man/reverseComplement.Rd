% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnaManipulate.R
\name{reverseComplement}
\alias{reverseComplement}
\title{Reverse-complementation of DNA}
\usage{
reverseComplement(nuc.sequences, reverse = TRUE)
}
\arguments{
\item{nuc.sequences}{Character vector containing the nucleotide sequences.}

\item{reverse}{Logical indicating if complement should be reversed.}
}
\value{
A character vector of reverse-complemented sequences.
}
\description{
The standard reverse-complement of nucleotide sequences.
}
\details{
With \samp{reverse = FALSE} the DNA sequence is only complemented, not reversed.

This function will handle the IUPAC ambiguity symbols, i.e. \samp{R} is
reverse-complemented to \samp{Y} etc.
}
\examples{
fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.ffn")
fa <- readFasta(fa.file)
reverseComplement(fa$Sequence)

#' # Or, make use of dplyr to manipulate tables
readFasta(fa.file) \%>\%
  mutate(RevComp = reverseComplement(Sequence)) -> fa.tbl

}
\seealso{
\code{\link{iupac2regex}}, \code{\link{regex2iupac}}.
}
\author{
Lars Snipen and Kristian Hovde Liland.
}

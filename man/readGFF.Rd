% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gff.R
\name{readGFF}
\alias{readGFF}
\alias{writeGFF}
\title{Reading and writing GFF-tables}
\usage{
readGFF(in.file)
writeGFF(gff.table, out.file)
}
\arguments{
\item{in.file}{Name of file with a GFF-table.}

\item{gff.table}{A table (\code{tibble}) with genomic features information.}

\item{out.file}{Name of file.}
}
\value{
\code{readGFF} returns a \code{gff.table} with the columns described above.

\code{writeGFF} writes the supplied \code{gff.table} to a text-file.
}
\description{
Reading or writing a GFF-table from/to file.
}
\details{
A GFF-table is simply a \code{\link{tibble}} with columns
adhering to the format specified by the GFF3 format, see
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details. There is
one row for each feature.

The following columns should always be in a full \code{gff.table} of the GFF3 format:
\itemize{
  \item Seqid. A unique identifier of the genomic sequence on which the feature resides.
  \item Source. A description of the procedure that generated the feature, e.g. \code{"R-package micropan::findOrfs"}.
  \item Type The type of feature, e.g. \code{"ORF"}, \code{"16S"} etc.
  \item Start. The leftmost coordinate. This is the start if the feature is on the Sense strand, but
  the end if it is on the Antisense strand.
  \item End. The rightmost coordinate. This is the end if the feature is on the Sense strand, but
  the start if it is on the Antisense strand.
  \item Score. A numeric score (E-value, P-value) from the \code{Source}. 
  \item Strand. A \code{"+"} indicates Sense strand, a \code{"-"} Antisense.
  \item Phase. Only relevant for coding genes. the values 0, 1 or 2 indicates the reading frame, i.e. 
  the number of bases to offset the \code{Start} in order to be in the reading frame.
  \item Attributes. A single string with semicolon-separated tokens prociding additional information.
}
Missing values are described by \code{"."} in the GFF3 format. This is also done here, except for the
numerical columns Start, End, Score and Phase. Here \code{NA} is used, but this is replaced by
\code{"."} when writing to file.

The \code{readGFF} function will also read files where sequences in FASTA format are added after
the GFF-table. This file section must always start with the line \code{##FASTA}. This fasta object
is added to the GFF-table as an attribute (use \code{attr(gff.tbl, "FASTA")} to retrieve it).
}
\examples{
# Using a GFF file in this package
gff.file <- file.path(path.package("microseq"),"extdata","small.gff")

# Reading gff-file
gff.tbl <- readGFF(gff.file)

}
\seealso{
\code{\link{findOrfs}}, \code{\link{lorfs}}.
}
\author{
Lars Snipen and Kristian Hovde Liland.
}

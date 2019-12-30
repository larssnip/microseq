#' @name gff2fasta
#' @title Retrieving annotated sequences
#' 
#' @description Retrieving from a genome the sequences specified in a \code{gff.table}.
#' 
#' @usage gff2fasta(gff.table, genome)
#' 
#' @param gff.table A \code{gff.table} (\code{tibble}) with genomic features information.
#' @param genome A fasta object (\code{tibble}) with the genome sequence(s).
#' 
#' @details Each row in \code{gff.table} (see \code{\link{readGFF}}) describes a genomic feature
#' in the \code{genome}, which is a \code{\link{tibble}} with columns \samp{Header} and
#' \samp{Sequence}. The information in the columns Seqid, Start, End and Strand are used to
#' retrieve the sequences from the \samp{Sequence} column of \code{genome}. Every Seqid in
#' the \code{gff.table} must match the first token in one of the \samp{Header} texts, in
#' order to retrieve from the correct \samp{Sequence}. 
#' 
#' @return A fasta object with one row for each row in \code{gff.table}. 
#' The \code{Header} for each sequence is a summary of the information in the
#' corresponding row of \code{gff.table}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{findOrfs}}.
#' 
#' @examples
#' # Using two files in this package
#' xpth <- file.path(path.package("microseq"),"extdata")
#' gff.file <- file.path(xpth,"small.gff")
#' genome.file <- file.path(xpth,"small_genome.fasta")
#' 
#' # Reading the genome first
#' genome <- readFasta(genome.file)
#' 
#' # Retrieving sequences
#' gff.table <- readGFF(gff.file)
#' fasta.tbl <- gff2fasta(gff.table, genome)
#' 
#' # Alternative, using piping
#' readGFF(gff.file) %>% gff2fasta(genome) -> fasta.tbl
#' 
#' @useDynLib microseq
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr mutate rename right_join select if_else %>% 
#' @importFrom stringr word str_sub str_c
#' 
#' @export gff2fasta
#' 
gff2fasta <- function(gff.table, genome){
  genome %>% 
    mutate(Header = word(Header, 1, 1)) %>% 
    rename(Seqid = Header, Gseq = Sequence) %>% 
    right_join(gff.table, by = "Seqid") %>% 
    mutate(Sequence = str_sub(Gseq, Start, End)) %>% 
    mutate(Sequence = if_else(Strand == "+", Sequence, reverseComplement(Sequence))) %>% 
    mutate(Header = str_c("Seqid=", Seqid, ";Start=", Start, ";End=", End, ";Strand=", Strand)) %>% 
    select(Header, Sequence) -> fsa
  class(fsa) <- c("Fasta", "data.frame")
  return(fsa)
}



#' @name readGFF
#' @title Reading and writing GFF-tables
#' @aliases readGFF writeGFF
#' 
#' @description Reading or writing a GFF-table from/to file.
#' 
#' @usage readGFF(in.file)
#' writeGFF(gff.table, out.file)
#' 
#' @param in.file Name of file with a GFF-table.
#' @param gff.table A table (\code{tibble}) with genomic features information.
#' @param out.file Name of file.
#' 
#' @details A GFF-table is simply a \code{\link{tibble}} with columns
#' adhering to the format specified by the GFF3 format, see
#' https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details. There is
#' one row for each feature.
#' 
#' The following columns should always be in a full \code{gff.table} of the GFF3 format:
#' \itemize{
#'   \item Seqid. A unique identifier of the genomic sequence on which the feature resides.
#'   \item Source. A description of the procedure that generated the feature, e.g. \code{"R-package micropan::findOrfs"}.
#'   \item Type The type of feature, e.g. \code{"ORF"}, \code{"16S"} etc.
#'   \item Start. The leftmost coordinate. This is the start if the feature is on the Sense strand, but
#'   the end if it is on the Antisense strand.
#'   \item End. The rightmost coordinate. This is the end if the feature is on the Sense strand, but
#'   the start if it is on the Antisense strand.
#'   \item Score. A numeric score (E-value, P-value) from the \code{Source}. 
#'   \item Strand. A \code{"+"} indicates Sense strand, a \code{"-"} Antisense.
#'   \item Phase. Only relevant for coding genes. the values 0, 1 or 2 indicates the reading frame, i.e. 
#'   the number of bases to offset the \code{Start} in order to be in the reading frame.
#'   \item Attributes. A single string with semicolon-separated tokens prociding additional information.
#' }
#' Missing values are described by \code{"."} in the GFF3 format. This is also done here, except for the
#' numerical columns Start, End, Score and Phase. Here \code{NA} is used, but this is replaced by
#' \code{"."} when writing to file.
#' 
#' The \code{readGFF} function will also read files where sequences in FASTA format are added after
#' the GFF-table. This file section must always start with the line \code{##FASTA}. This fasta object
#' is added to the GFF-table as an attribute (use \code{attr(gff.tbl, "FASTA")} to retrieve it).
#' 
#' @return \code{readGFF} returns a \code{gff.table} with the columns described above.
#' 
#' \code{writeGFF} writes the supplied \code{gff.table} to a text-file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{findOrfs}}, \code{\link{lorfs}}.
#' 
#' @examples
#' # Using a GFF file in this package
#' xpth <- file.path(path.package("microseq"),"extdata")
#' gff.file <- file.path(xpth,"small.gff")
#' 
#' # Reading gff-file
#' gff.tbl <- readGFF(gff.file)
#' 
#' @importFrom stringr str_split str_c
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate_at %>% 
#' 
#' @export readGFF writeGFF
#' 
readGFF <- function(in.file){
  fil <- file(in.file, open = "rt")
  lines <- readLines(fil)
  close(fil)
  fasta.idx <- grep("##FASTA", lines)
  cn <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
  if(length(lines) > 1){
    if(length(fasta.idx) > 0){
      lns1 <- lines[1:fasta.idx]
      gff.table <- as_tibble(str_split(lns1[!grepl("^#", lns1)], pattern = "\t", simplify = T))
      if(ncol(gff.table) != 9 ) stop("Table must have 9 tab-separated columns, this one has", ncol(gff.table))
      colnames(gff.table) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
      gff.table %>% 
        mutate_at(c("Start", "End", "Score", "Phase"), as.numeric) -> gff.table
      lns2 <- lines[(fasta.idx+1):length(lines)]
      idx <- c(grep("^>", lns2), length(lns2) + 1)
      fsa <- tibble(Header = gsub("^>", "", lns2[idx[1:(length(idx)-1)]]),
                    Sequence = sapply(1:(length(idx)-1), function(ii){
                      str_c(lns2[(idx[ii]+1):(idx[ii+1]-1)], collapse = "")
                    }))
      attr(gff.table, "FASTA") <- fsa
    } else {
      gff.table <- as_tibble(str_split(lines[!grepl("^#", lines)], pattern = "\t", simplify = T))
      if(ncol(gff.table) != 9 ) stop("Table must have 9 tab-separated columns, this one has", ncol(gff.table))
      colnames(gff.table) <- cn
      gff.table %>% 
        mutate_at(c("Start", "End", "Score", "Phase"), as.numeric) -> gff.table
    }
  } else {
    gff.table <- tibble("Seqid" = character(0),
                        "Source" = character(0),
                        "Type" = character(0),
                        "Start" = numeric(0),
                        "End" = numeric(0),
                        "Score" = numeric(0),
                        "Strand" = character(0),
                        "Phase" = numeric(0),
                        "Attributes" = character(0))
  }
  return(gff.table)
}
writeGFF <- function(gff.table, out.file){
  line1 <- c("##gff-version 3")
  sapply(1:nrow(gff.table), function(i){str_c(gff.table[i,], collapse = "\t")}) %>% 
    str_replace_all("\tNA\t", "\t.\t") %>% 
    str_replace_all("\tNA$", "\t.") -> lines
  writeLines(c(line1, lines), con = out.file)
  return(NULL)
}

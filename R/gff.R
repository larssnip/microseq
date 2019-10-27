#' @name gff2fasta
#' @title Retrieving annotated sequences
#' 
#' @description Retrieving from a genome the sequences specified in a \code{gff.table}.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' @param genome A \code{\link{Fasta}} object with the genome sequence(s).
#' 
#' @details Each row in \code{gff.table} (see \code{\link{readGFF}}) describes a genomic feature
#' in the \code{genome}. The information in the columns Seqid, Start, End and Strand are used to retrieve
#' the sequences from \code{genome$Sequence}. Every Seqid in the \code{gff.table}
#' must match the first token in one of the \code{genome$Header} texts.
#' 
#' @return A \code{\link{Fasta}} object with one row for each row in \code{gff.table}. 
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
#' # Reading
#' gff.table <- readGFF(gff.file)
#' genome <- readFasta(genome.file)
#' 
#' # Retrieving sequences
#' fasta.obj <- gff2fasta(gff.table,genome)
#' summary(fasta.obj)
#' plot(fasta.obj)
#' 
#' @useDynLib microseq
#' @importFrom Rcpp evalCpp
#' 
#' @export gff2fasta
#' 
gff2fasta <- function(gff.table, genome){
  utag <- unique(gff.table$Seqid)
  tagz <- unique(sapply(strsplit(genome$Header, split = " "), function(x){x[1]}))
  if(sum(is.na(match(utag,tagz))) > 0)
    stop( "Seqid's in the gff.table do not match the genome$Header" )
  
  n <- nrow(gff.table)
  fobj <- data.frame(Header=gffSignature(gff.table), 
                     Sequence=rep("", n),
                     stringsAsFactors = F)
  strnd <- rep(1, n)
  strnd[gff.table$Strand == "-"] <- -1
  gff.table$Strand <- strnd
  for(i in 1:length(tagz)){
    idx <- which(gff.table$Seqid == tagz[i])
    seq <- extractSeq(genome$Sequence[i], 
                      gff.table$Start[idx], 
                      gff.table$End[idx], 
                      gff.table$Strand[idx])
    idd <- which(gff.table$Strand[idx] < 0)
    seq[idd] <- reverseComplement(seq[idd])
    fobj$Sequence[idx] <- seq
  }
  class(fobj) <- c("Fasta", "data.frame")
  return(fobj)
}



#' @name gffSignature
#' @title GFF signature text
#' 
#' @description Making a signature text from \code{gff.table} data.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' 
#' @details For each row in \code{gff.table} a text is created by pasting these
#' data together, adding some explanatory text. This function is used by \code{link{gff2fasta}}
#' to create the Header-lines for the \code{\link{Fasta}} object when retrieving the sequences.
#' 
#' @return A vector of texts, one for each row in \code{gff.table}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{findOrfs}}, \code{\link{gff2fasta}}.
#' 
#' @examples # See the example in the Help-file for readGFF.
#' 
#' @useDynLib microseq
#' 
#' @export gffSignature
#' 
gffSignature <- function(gff.table){
  desc <- paste0("Seqid=", gff.table$Seqid,
                 ";Source=", gff.table$Source,
                 ";Type=", gff.table$Type,
                 ";Start=", gff.table$Start,
                 ";End=", gff.table$End,
                 ";Score=", gff.table$Score,
                 ";Strand=", gff.table$Strand,
                 ";Phase=", gff.table$Phase,
                 ";Attributes=", gff.table$Attributes)
  return(desc)
}




#' @name readGFF
#' @title Reading and writing GFF-tables
#' @aliases readGFF writeGFF
#' 
#' @description Reading or writing a GFF-table from/to file.
#' 
#' @usage readGFF(in.file)
#' writeGFF(gff.tbl, out.file)
#' 
#' @param in.file Name of file with a GFF-table.
#' @param gff.tbl A table (\code{data.frame} or \code{tibble}) with genomic features information.
#' @param out.file Name of file.
#' 
#' @details A GFF-table is simply a \code{\link{data.frame}} or \code{\link{tibble}} with columns
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
#' The \code{readGFF} function will also read files where sequences in FASTA format are added after the GFF-table.
#' This file section must always start with the line \code{##FASTA}. This \code{\link{Fasta}} object is added to
#' the GFF-table as an attribute (use \code{attr(gff.tbl, "Fasta")} to retrieve it).
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
#' # Reading gff-file and print signature
#' gff.tbl <- readGFF(gff.file)
#' print(gffSignature(gff.tbl))
#' 
#' @importFrom stringr str_split str_c
#' @importFrom tibble tibble as_tibble
#' 
#' @export readGFF writeGFF
#' 
readGFF <- function(in.file){
  fil <- file(in.file, open = "rt")
  lines <- readLines(fil)
  close(fil)
  fasta.idx <- grep("##FASTA", lines)
  if(length(lines) > 1){
    if(length(fasta.idx) > 0){
      lns1 <- lines[1:fasta.idx]
      tbl <- as_tibble(str_split(lns1[!grepl("^#", lns1)], pattern = "\t", simplify = T))
      if(ncol(tbl) != 9 ) stop("Table must have 9 tab-separated columns, this one has", ncol(tbl))
      colnames(tbl) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
      w <- options()$warn
      options(warn = -1)
      tbl$Start <- as.numeric(tbl$Start)
      tbl$End   <- as.numeric(tbl$End)
      tbl$Score <- as.numeric(tbl$Score)
      tbl$Phase <- as.numeric(tbl$Phase)
      options(warn = w)
      
      lns2 <- lines[(fasta.idx+1):length(lines)]
      idx <- c(grep("^>", lns2), length(lns2) + 1)
      fsa <- tibble(Header = gsub("^>", "", lns2[idx[1:(length(idx)-1)]]),
                    Sequence = sapply(1:(length(idx)-1), function(ii){
                      str_c(lns2[(idx[ii]+1):(idx[ii+1]-1)], collapse = "")
                    }))
      attr(tbl, "Fasta") <- fsa
    } else {
      tbl <- as_tibble(str_split(lines[!grepl("^#", lines)], pattern = "\t", simplify = T))
      if(ncol(tbl) != 9 ) stop("Table must have 9 tab-separated columns, this one has", ncol(tbl))
      colnames(tbl) <- c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
      w <- options()$warn
      options(warn = -1)
      tbl$Start <- as.numeric(tbl$Start)
      tbl$End   <- as.numeric(tbl$End)
      tbl$Score <- as.numeric(tbl$Score)
      tbl$Phase <- as.numeric(tbl$Phase)
      options(warn = w)
    }
  } else {
    tbl <- tibble("Seqid" = NULL,
                  "Source" = NULL,
                  "Type" = NULL,
                  "Start" = NULL,
                  "End" = NULL,
                  "Score" = NULL,
                  "Strand" = NULL,
                  "Phase" = NULL,
                  "Attributes" = NULL)
  }
  return(tbl)
}
writeGFF <- function(gff.tbl, out.file){
  line1 <- c("##gff-version 3")
  lines <- sapply(1:nrow(gff.tbl), function(i){paste(gff.tbl[i,], collapse = "\t")})
  lines <- gsub("\tNA\t", "\t.\t", lines)
  lines <- gsub("\tNA$", "\t.", lines)
  writeLines(c(line1, lines), con = out.file)
  return(paste("gff.table written to", out.file))
}

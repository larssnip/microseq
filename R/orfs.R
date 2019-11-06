#' @name findOrfs
#' @title Finding ORFs in genomes
#' 
#' @description Finds all ORFs in prokaryotic genome sequences.
#' 
#' @param genome A \code{\link{Fasta}} object with the genome sequence(s).
#' @param circular Logical indicating if the genome sequences are completed, circular sequences.
#' @param trans.tab Translation table
#' 
#' @details A prokaryotic Open Reading Frame (ORF) is defined as a subsequence starting with a  start-codon
#' (ATG, GTG or TTG), followed by an integer number of triplets (codons), and ending with a stop-codon (TAA,
#' TGA or TAG, unless \code{trans.tab = 4}, see below). This function will locate all such ORFs in a genome.
#' 
#' The argument \code{genome} will typically have several sequences (chromosomes/plasmids/scaffolds/contigs).
#' It is vital that the \emph{first token} (characters before first space) of every \code{genome$Header} is
#' unique, since this will be used to identify these genome sequences in the output.
#' 
#' An alternative translation table may be specified, and as of now the only alternative implemented is table 4.
#' This means codon TGA is no longer a stop, but codes for Tryptophan. This coding is used by some bacteria
#' (e.g. Mycoplasma, Mesoplasma).
#' 
#' Note that for any given stop-codon there are usually multiple start-codons in the same reading
#' frame. This function will return all, i.e. the same stop position may appear multiple times. If
#' you want ORFs with the most upstream start-codon only (LORFs), then filter the output from this function
#' with \code{\link{lorfs}}.
#' 
#' By default the genome sequences are assumed to be linear, i.e. contigs or other incomplete fragments
#' of a genome. In such cases there will usually be some truncated ORFs at each end, i.e. ORFs where either
#' the start- or the stop-codon is lacking. In the \code{orf.table} returned by this function this is marked in the
#' Attributes column. The texts "Truncated=10" or "Truncated=01" indicates truncated at 
#' the beginning or end of the genomic sequence, respectively. If the supplied \code{genome} is a completed genome, with 
#' circular chromosome/plasmids, set the flag \code{circular = TRUE} and no truncated ORFs will be listed.
#' In cases where an ORF runs across the origin of a circular genome sequences, the Stop coordinate will be
#' larger than the length of the genome sequence. This is in line with the specifications of the GFF3 format, where 
#' a Start cannot be larger than the corresponding End.
#' 
#' @return This function returns an \code{orf.table}, which is simply a \code{\link{tibble}} with columns
#' adhering to the GFF3 format specifications (a \code{gff.table}), see \code{\link{readGFF}}. If you want to retrieve
#' the ORF sequences, use \code{\link{gff2fasta}}.
#' 
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}, \code{\link{lorfs}}.
#' 
#' @examples
#' # Using a genome file in this package
#' xpth <- file.path(path.package("microseq"),"extdata")
#' genome.file <- file.path(xpth,"small_genome.fasta")
#' 
#' # Reading genome and finding orfs
#' genome <- readFasta(genome.file)
#' orf.tbl <- findOrfs(genome)
#' 
#' # Pipeline for finding LORFs of minimum length 100 amino acids
#' # and collecting their sequences from the genome
#' findOrfs(genome) %>% 
#'  lorfs() %>% 
#'  filter(orfLength(., aa = T) > 100) %>% 
#'  gff2fasta(genome) -> lorf.fasta
#' 
#' @useDynLib microseq
#' @importFrom Rcpp evalCpp
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate right_join filter bind_rows select if_else
#' @importFrom stringr word str_length
#' 
#' @export findOrfs
#' 
findOrfs <- function(genome, circular = F, trans.tab = 11){
  genome %>% 
    mutate(Header = word(Header, 1, 1)) %>% 
    mutate(Length = str_length(Sequence)) -> genome
  if(length(unique(genome$Header)) != length(genome$Header)) stop("First token in the Headers must be unique!")
  ORF_index(genome$Header, genome$Sequence, trans.tab) %>% 
    as_tibble() %>% 
    mutate(Seqid = as.character(Seqid)) -> orf.table
  if(circular){
    orf.table %>% 
      right_join(genome, by = c("Seqid" = "Header")) -> ott
    otn <- circularize(ott, trans.tab)
    orf.table %>% 
      filter(Truncated == 0) %>% 
      bind_rows(otn) -> orf.table
  }
  orf.table %>% 
    mutate(Strand = if_else(Strand > 0, "+", "-")) %>% 
    mutate(Attributes = "Truncated=00") %>% 
    mutate(Attributes = if_else(Truncated > 0, "Truncated=10", Attributes)) %>% 
    mutate(Attributes = if_else(Truncated < 0, "Truncated=01", Attributes)) %>% 
    mutate(Source = "micropan::findOrfs") %>% 
    mutate(Type = "ORF", Score = NA, Phase = 0) %>% 
    select(Seqid, Source, Type, Start, End, Score, Strand, Phase, Attributes) -> orf.table
  return(orf.table)
}

# Local function
circularize <- function(ot, trans.tab){
  ugs <- unique(ot$Seqid)
  otn <- NULL
  for(i in 1:length(ugs)){
    ot %>% 
      filter(Seqid == ugs[i]) -> otg
    if(max(otg$Truncated) != min(otg$Truncated)){
      gseq.pre  <- str_sub(otg$Sequence[1], 1, 10000)
      gseq.post <- str_sub(otg$Sequence[1], -10000, -1)
      dd <- otg$Length[1] - str_length(gseq.post)
      ORF_index(otg$Seqid[1], str_c(gseq.post, gseq.pre), trans.tab) %>% 
        as_tibble() %>% 
        mutate(Seqid = as.character(Seqid)) %>% 
        filter(Truncated == 0) %>% 
        mutate(Start = Start + dd) %>% 
        mutate(End = End + dd) %>% 
        filter(Start < otg$Length[1] & End > otg$Length[1]) %>% 
        bind_rows(otn) -> otn
    }
  }
  return(otn)
}




#' @name orfLength
#' @title Length of ORF
#' 
#' @description Computing the lengths of all ORFs in an \code{orf.table}.
#' 
#' @param orf.table A \code{\link{gff.table}} with ORF information.
#' @param aa Logical, length in amino acids instead of bases.
#' 
#' @details Computes the length of an ORF in bases, including the stop codon. However, if \code{aa = TRUE},
#' then the length is in amino acids after translation. This aa-length is the base-length divided by 3 and minus 1,
#' unless the ORF is truncated and lacks a stop codon.
#' 
#' @return A vector of integers.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{findOrfs}}.
#' 
#' @examples # See the example in the Help-file for findOrfs.
#' 
#' @importFrom dplyr mutate if_else %>% 
#' 
#' @export orfLength
#' 
orfLength <- function(orf.table, aa = FALSE){
  orf.table %>% 
    mutate(Length = abs(Start - End) + 1) -> orf.table
  if(aa){
    orf.table %>% 
      mutate(Length = Length / 3) %>% 
      mutate(Length = if_else((str_detect(Attributes, "Truncated=01") & Strand == "+")
                              |(str_detect(Attributes, "Truncated=10") & Strand == "-"),
                              Length, Length - 1)) -> orf.table
  }
  return(orf.table$Length)
}


#' @name lorfs
#' @title Longest ORF
#' 
#' @description Filtering an \code{orf.table} with ORF information to keep only the LORFs.
#' 
#' @param orf.table A \code{tibble} with the nine columns of the GFF-format (see \code{\link{findOrfs}}).
#' 
#' @details For every stop-codon there are usually multiple possible start-codons in the same reading
#' frame (nested ORFs). The LORF (Longest ORF) is defined as the longest of these nested ORFs,
#' i.e. the ORF starting at the most upstream start-codon matching the stop-codon.
#' 
#' @return A \code{\link{tibble}} with a subset of the rows of the argument \code{orf.table}. 
#' After this filtering the Type variable in \code{orf.table} is changed to \code{"LORF"}. If you want to
#' retrieve the LORF sequences, use \code{\link{gff2fasta}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{findOrfs}}, \code{\link{gff2fasta}}.
#' 
#' @examples # See the example in the Help-file for findOrfs.
#' 
#' @importFrom dplyr mutate arrange distinct select %>% 
#' 
#' @export lorfs
#' 
lorfs <- function(orf.table){
  orf.table %>% 
    mutate(Length = orfLength(.)) %>% 
    arrange(desc(Length)) %>% 
    mutate(End.end = if_else(Strand == "+", End, Start)) %>% 
    mutate(Signature = str_c(Seqid, End.end, Strand)) %>% 
    distinct(Signature, .keep_all = T) %>% 
    mutate(Type = "LORF") %>% 
    select(-End.end, -Signature, -Length) -> lorf.tbl
  return(lorf.tbl)
}

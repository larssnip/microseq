#' @name msalign
#' @title Multiple alignment
#' 
#' @description Quickly computing a smallish multiple sequence alignment.
#' 
#' @param fsa.tbl A fasta object (data.frame or tibble) with input sequences.
#' @param machine Function that does the 'dirty work'.
#' 
#' @details This function computes a multiple sequence alignment given a set of sequences in
#' a fasta object, see \code{\link{readFasta}} for more on fasta objects.
#' 
#' It is merely a wrapper for the function named in \code{machine} to avoid explicit writing
#' and reading of files. This function should only be used for small data sets, since no result
#' files are stored. For heavier jobs, use the \code{machine} function directly.
#' 
#' At present, the only \code{machine} function implemented is \code{\link{muscle}}, but other
#' third-party \code{machine}s may be included later.
#' 
#' Note that this function will run \code{\link{muscle}} with default settings, which is fine
#' for small data sets.
#' 
#' @return Results are returned as a fasta object, i.e. a tibble with columns
#' \samp{Header} and \samp{Sequence}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{muscle}}, \code{\link{msaTrim}}.
#' 
#' @examples
#' \dontrun{
#' prot.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.faa")
#' faa <- readFasta(prot.file)
#' msa <- msalign(faa)
#' }
#' 
#' @importFrom stats rnorm
#' @importFrom stringr str_c
#' 
#' @export msalign
#' 
msalign <- function(fsa.tbl, machine = "microseq::muscle"){
  available.external(machine)
  in.fil <- file.path(getwd(), "in_file.fa")
  ut.fil <- file.path(getwd(), "out_file.fa")
  writeFasta(fsa.tbl, out.file = in.fil)
  cmd <- str_c(machine, "('", in.fil, "','", ut.fil, "', quiet = T)")
  eval(parse(text = cmd))
  msa.tbl <- readFasta(ut.fil)
  ok <- file.remove(c(in.fil, ut.fil))
  return(msa.tbl)
}


#' @name msaTrim
#' @title Trimming multiple sequence alignments
#' 
#' @description Trimming a multiple sequence alignment by discarding columns with too many gaps.
#' 
#' @param msa A fasta object containing a multiple alignment.
#' @param gap.end Fraction of gaps tolerated at the ends of the alignment (0-1).
#' @param gap.mid Fraction of gaps tolerated inside the alignment (0-1).
#' 
#' @details A multiple alignment is trimmed by removing columns with too many indels (gap-symbols). Any 
#' columns containing a fraction of gaps larger than \code{gap.mid} are discarded. For this reason, \code{gap.mid}
#' should always be farily close to 1.0 therwise too many columns may be discarded, destroying the alignment.
#' 
#' Due to the heuristics of multiple alignment methods, both ends of the alignment tend to be uncertain and most
#' of the trimming should be done at the ends. Starting at each end, columns are discarded as long as their fraction of gaps
#' surpasses \code{gap.end}. Typically \code{gap.end} can be much smaller than \code{gap.mid}, but if 
#' set too low you risk that all columns are discarded!
#' 
#' @return The trimmed alignment is returned as a fasta object.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{muscle}}, \code{\link{msalign}}.
#' 
#' @examples 
#' msa.file <- file.path(path.package("microseq"),"extdata", "small.msa")
#' msa <- readFasta(msa.file)
#' print(str_length(msa$Sequence))
#' msa.trimmed <- msaTrim(msa)
#' print(str_length(msa.trimmed$Sequence))
#' msa.mat <- msa2mat(msa)  # for use with ape::as.DNAbin(msa.mat)
#' 
#' @importFrom stringr str_length str_split
#' 
#' @export msaTrim
#' 
msaTrim <- function(msa, gap.end = 0.5, gap.mid = 0.9){
  nc <- unique(str_length(msa$Sequence))
  if(length(nc) != 1) stop("This is not a multiple alignment, sequences have different lengths!")
  cmat <- str_split(msa$Sequence, pattern = "", simplify = T)
  gap.frac <- colSums(cmat == "-")/nrow(cmat)
  if(min(gap.frac) > min(gap.end, gap.mid)){
    warning("All positions have more gaps than arguments allow, please increase gap.end value")
    msa$Sequence <- ""
    return(msa)
  }
  idx.keep <- (gap.frac <= gap.mid)
  j <- 1
  while(gap.frac[j] > gap.end){
    idx.keep[j] <- F
    j <- j + 1
  }
  j <- nc
  while(gap.frac[j] > gap.end){
    idx.keep[j] <- F
    j <- j - 1
  }
  msa$Sequence <- apply(cmat[,idx.keep], 1, paste, collapse="")
  return(msa)
}



#' @name msa2mat
#' @title Convert alignment to matrix
#' 
#' @description Converts a FASTA formatted multiple alignment to a matrix.
#' 
#' @param msa A fasta object with a multiple alignment, see \code{\link{msalign}}`.
#' 
#' @details This function converts the fasta object \code{msa}, containing a multiple alignment,
#' to a matrix. This means each position in the alignment is a column in the matrix, and the
#' content of the \samp{Header} column of \code{msa} is used as rownames of theh matrix.
#' 
#' Such a matrix is useful for conversion to a \code{DNAbin} object that is used by the \code{ape}
#' package for re-constructing phylogenetic trees. 
#' 
#' @return A \code{matrix} where each row is a vector of aligned bases/amino acids.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{msalign}}, \code{\link{readFasta}}.
#' 
#' @examples
#' msa.file <- file.path(path.package("microseq"),"extdata", "small.msa")
#' msa <- readFasta(msa.file)
#' msa.mat <- msa2mat(msa)  # to use with ape::as.DNAbin(msa.mat)
#' 
#' @importFrom stringr str_split
#' 
#' @export msa2mat
#' 
msa2mat <- function(msa){
  msa.mat <- str_split(msa$Sequence, pattern = "", simplify = T)
  rownames(msa.mat) <- msa$Header
  return(msa.mat)
}




#' @name backTranslate
#' @title Replace amino acids with codons
#' 
#' @description Replaces aligned amino acids with their original codon triplets.
#' 
#' @param aa.msa A fasta object with a multiple alignment, see \code{\link{msalign}}.
#' @param nuc.ffn A fasta object with the coding sequences, see \code{\link{readFasta}}.
#' 
#' @details This function replaces the aligned amino acids in \code{aa.msa} with their
#' original codon triplets. This is possible only when the nucleotide sequences in \code{nuc.ffn}
#' are the exact nucleotide sequences behind the protein sequences that are aligned in \code{aa.msa}.
#' 
#' It is required that the first token of the \samp{Header} lines is identical for a protein sequence
#' in \code{aa.msa} and its nucleotide version in \samp{nuc.ffn}, otherwise it is impossible to
#' match them. Thus, they may not appear in the same order in the two input fasta objects.
#' 
#' When aligning coding sequences, one should in general always align their protein sequences, to
#' keep the codon structure, and then use \code{\link{backTranslate}} to convert this into a 
#' nucleotide alignment, if required.
#' 
#' If the nuclotide sequences contain the stop codons, these will be removed.
#' 
#' @return A fasta object similar to \code{aa.msa}, but where each amino acid has been replace by 
#' its corresponding codon. All gaps \samp{"-"} are replaced by triplets \samp{"---"}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{msalign}}, \code{\link{readFasta}}.
#' 
#' @examples
#' msa.file <- file.path(file.path(path.package("microseq"),"extdata"), "small.msa")
#' aa.msa <- readFasta(msa.file)
#' nuc.file <- file.path(file.path(path.package("microseq"),"extdata"), "small.ffn")
#' nuc <- readFasta(nuc.file)
#' nuc.msa <- backTranslate(aa.msa, nuc)
#' 
#' @importFrom stringr word str_remove str_sub str_length str_c
#' @importFrom dplyr slice mutate %>% 
#' @importFrom rlang .data
#' 
#' @export backTranslate
#' 
backTranslate <- function(aa.msa, nuc.ffn){
  ixx <- match(word(aa.msa$Header, 1, 1), word(nuc.ffn$Header, 1, 1))
  if(sum(is.na(ixx)) > 0) stop("The Headers of aa.msa and nuc.ffn must have matching first tokens")
  nuc.ffn %>% 
    slice(ixx) %>% 
    mutate(Sequence = str_remove(.data$Sequence, "TGA$|TAG$|TAA$")) -> nuc.ffn
  aa.msa %>% 
    select(.data$Header, .data$Sequence) -> nuc.msa
  M.prot <- msa2mat(aa.msa)
  for(i in 1:nrow(aa.msa)){
    codons <- str_sub(nuc.ffn$Sequence[i],
                      start = seq(1, str_length(nuc.ffn$Sequence[i]), 3),
                      end = seq(3, str_length(nuc.ffn$Sequence[i]), 3))
    idx <- which(M.prot[i,] != "-")
    if(length(codons) != length(idx)) stop("Number of amino acids in sequence ", i, " does not match the number of codons in corresponding nucleotide sequence")
    M.nuc <- rep("---", ncol(M.prot))
    M.nuc[idx] <- codons
    nuc.msa$Sequence[i] <- str_c(M.nuc, collapse = "")
  }
  return(nuc.msa)
}



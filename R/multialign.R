#' @name msalign
#' @title Multiple alignment
#' 
#' @description Quickly computing a smallish multiple sequence alignment.
#' 
#' @usage msalign(fdta, machine = "muscle")
#' 
#' @param fdta A fasta object (\code{tibble}) with input sequences.
#' @param machine Function that does the 'dirty work'.
#' 
#' @details This function computes a multiple sequence alignment given a set of sequences in a fasta object,
#' see \code{\link{readFasta}} for more on fasta objects.
#' 
#' It is merely a wrapper for the function named in \code{machine} to avoid explicit writing and reading of files.
#' This function should only be used for small data sets, since no result files are stored. For heavier jobs,
#' use the \code{machine} function directly.
#' 
#' At present, the only \code{machine} function implemented is \code{\link{muscle}}, but other third-party \code{machine}s 
#' may be included later.
#' 
#' Note that this function will run \code{\link{muscle}} with default settings, which is fine for small data sets.
#' 
#' @return Results are returned as a fasta object, i.e. a \code{\link{tibble}} with columns \samp{Header} and \samp{Sequence}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{muscle}}, \code{\link{msaTrim}}.
#' 
#' @importFrom stats rnorm
#' 
#' @examples
#' \dontrun{
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(fa.file)
#' msa <- msalign(fdta)
#' }
#' 
#' @export msalign
#' 
msalign <- function(fdta, machine = "muscle"){
  available.external(machine)
  num <- substr(as.character(rnorm(1)), 4, 7)
  in.fil <- file.path(getwd(), paste0("in", num, ".fasta"))
  ut.fil <- file.path(getwd(), paste0("ut", num, ".fasta", sep=""))
  writeFasta(fdta, out.file = in.fil)
  cmd <- paste0(machine, "('", in.fil, "','", ut.fil, "',quiet = T)")
  eval(parse(text = cmd))
  msa <- readFasta(ut.fil)
  file.remove(c(in.fil, ut.fil))
  return(msa)
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
#' @importFrom stringr str_length str_split
#' 
#' @examples 
#' \dontrun{
#' msa.file <- file.path(file.path(path.package("microseq"),"extdata"), "msa.fasta")
#' msa <- readFasta(msa.file)
#' print(str_length(msa$Sequence))
#' msa.trimmed <- msaTrim(msa)
#' print(str_length(msa.trimmed$Sequence))
#' msa.mat <- msa2mat(msa)  # for use with as.DNAbin(msa.mat)
#' }
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
#' @usage msa2mat(msa)
#' 
#' @param msa A fasta object with a multiple alignment, see \code{\link{msalign}}`.
#' 
#' @details This function converts the fasta object \code{msa}, containing a multiple alignment, to a matrix. This means
#' each position in the alignment is a column in the matrix, and the content of the \samp{Header} column of \code{msa}
#' is used as rownames of theh matrix.
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
#' \dontrun{
#' msa.file <- file.path(file.path(path.package("microseq"),"extdata"), "msa.fasta")
#' msa <- readFasta(msa.file)
#' msa.mat <- msa2mat(msa)
#' }
#' 
#' @keywords sequence FASTA
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

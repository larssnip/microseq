#' @name translate
#' @title Translation according to the standard genetic code
#' 
#' @description The translation from DNA(RNA) to amino acid sequence according to the standard genetic code.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' @param M.start A logical indicating if the amino acid sequence should start with M regardless of start codon.
#' @param no.stop A logical indicating if terminal stops (*) should be eliminated from the translated sequence
#' @param trans.tab Translation table, either 11 or 4
#' 
#' @details Codons are by default translated according to translation table 11, i.e. the possible start codons
#' are ATG, GTG or TTG and stop codons are TAA, TGA and TAG.  The only alternative implemented here is
#' translation table 4, which is used by some bacteria (e.g. Mycoplasma, Mesoplasma). If \code{trans.tab} is 4,
#' the stop codon TGA istranslated to W (Tryptophan).
#' 
#' @return A character vector of translated sequences.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.ffn")
#' fa <- readFasta(fa.file)
#' translate(fa$Sequence)
#' 
#' # Or, make use of dplyr to manipulate tables
#' readFasta(fa.file) %>%
#'   mutate(Protein = translate(Sequence)) -> fa.tbl
#' 
#' @importFrom stringr str_replace_all
#' 
#' @export translate
#' 
translate <- function(nuc.sequences, M.start = TRUE, no.stop = TRUE, trans.tab = 11){
  nuc.sequences <- str_replace_all(toupper(nuc.sequences), "U", "T")
  if(M.start) nuc.sequences <- str_replace_all(nuc.sequences, "^.{3}", "ATG")
  prot.sequences <- transl(nuc.sequences, trans.tab)
  if(no.stop) prot.sequences <- gsub("\\*$", "", prot.sequences)
  return(prot.sequences)
}



#' @name ctranslate
#' @title Translating codons into a single letter alphabet
#' 
#' @description Codons of a coding genes is transpalted into a single letter alphabet.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' 
#' @details Codons are translated into an alphabet of 64 single letters, the 26 lowercase
#' english letters, the 26 uppercase english letters, digits 1-9 as well as the symbols
#' "!", ":" and "?" for the three stop-codons. Thus, the translated sequence has the same 
#' length as the amino acid sequence (see \code{\link{microseq::translate}}) but the
#' alphabet is different.
#' 
#' Any letter in the DNA sequence that is not A, C, G or T will result in the corresponding 
#' codon symbol "*", indicating a non-valid codon.
#' 
#' @return A character vector of c-translated sequences.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.ffn")
#' fa <- readFasta(fa.file)
#' ctranslate(fa$Sequence)
#' 
#' # Or, make use of dplyr to manipulate tables
#' fa.tbl <- readFasta(fa.file) %>%
#'   mutate(Codon_seq = ctranslate(Sequence))
#' 
#' @importFrom stringr str_replace_all
#' 
#' @export ctranslate
#' 
ctranslate <- function(nuc.sequences){
  nuc.sequences <- str_replace_all(toupper(nuc.sequences), "U", "T")
  cds.sequences <- ctransl(nuc.sequences)
  return(cds.sequences)
}



#' @name reverseComplement
#' @title Reverse-complementation of DNA
#' 
#' @description The standard reverse-complement of nucleotide sequences.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' @param reverse Logical indicating if complement should be reversed.
#' 
#' @details With \samp{reverse = FALSE} the DNA sequence is only complemented, not reversed.
#' 
#' This function will handle the IUPAC ambiguity symbols, i.e. \samp{R} is
#' reverse-complemented to \samp{Y} etc.
#' 
#' @return A character vector of reverse-complemented sequences.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{iupac2regex}}, \code{\link{regex2iupac}}.
#' 
#' @examples 
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.ffn")
#' fa <- readFasta(fa.file)
#' reverseComplement(fa$Sequence)
#' 
#' #' # Or, make use of dplyr to manipulate tables
#' readFasta(fa.file) %>%
#'   mutate(RevComp = reverseComplement(Sequence)) -> fa.tbl
#' 
#' @export reverseComplement
#' 
reverseComplement <- function(nuc.sequences, reverse = TRUE){
  return(revComp(nuc.sequences, reverse))
}



#' @name iupac2regex
#' @title Ambiguity symbol conversion
#' @aliases iupac2regex regex2iupac
#' 
#' @description Converting DNA ambiguity symbols to regular expressions, and vice versa.
#' 
#' @usage iupac2regex(sequence)
#' regex2iupac(sequence)
#' 
#' @param sequence Character vector containing DNA sequences.
#' 
#' @details The DNA alphabet may contain ambiguity symbols, e.g. a W means either A or T.
#' When using a regular expression search, these letters must be replaced by the proper
#' regular expression, e.g. W is replaced by [AT] in the string. The \code{iupac2regex} makes this
#' translation, while \code{regex2iupac} converts the other way again (replace [AT] with W).
#' 
#' @return A string where the ambiguity symbol has been replaced by a regular expression
#' (\code{iupac2regex}) or a regular expression has been replaced by an ambiguity symbol
#' (\code{regex2iupac}).
#' 
#' @author Lars Snipen.
#' 
#' @examples
#' iupac2regex("ACWGT")
#' regex2iupac("AC[AG]GT")
#' 
#' @importFrom stringr str_replace_all fixed
#' 
#' @export iupac2regex
#' @export regex2iupac
#' 
iupac2regex <- function(sequence){
  IUPAC <- matrix(c("W","[AT]",
                    "S","[CG]",
                    "M","[AC]",
                    "K","[GT]",
                    "R","[AG]",
                    "Y","[CT]",
                    "B","[CGT]",
                    "D","[AGT]",
                    "H","[ACT]",
                    "V","[ACG]",
                    "N","[ACGT]"), ncol = 2, byrow = T)
  s <- str_replace_all(toupper(sequence), "X", "N")
  for(i in 1:nrow(IUPAC)) s <- str_replace_all(s, fixed(IUPAC[i,1]), IUPAC[i,2])
  return(s)
}
regex2iupac <- function(sequence){
  IUPAC <- matrix(c("W","[AT]",
                    "S","[CG]",
                    "M","[AC]",
                    "K","[GT]",
                    "R","[AG]",
                    "Y","[CT]",
                    "B","[CGT]",
                    "D","[AGT]",
                    "H","[ACT]",
                    "V","[ACG]",
                    "N","[ACGT]"), ncol = 2, byrow = T)
  s <- toupper(sequence)
  for(i in 1:nrow(IUPAC)) s <- str_replace_all(s, fixed(IUPAC[i,2]), IUPAC[i,1])
  return(s)
}


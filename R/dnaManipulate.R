#' @name translate
#' @title Translation according to the standard genetic code
#' 
#' @description The translation from DNA(RNA) to amino acid sequence according to the standard genetic code.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' @param M.start A logical indicating if the amino acid sequence should start with M regardless of 
#' start codon (ATG, GTG or TTG).
#' 
#' @details This function uses the Biostrings::translate function.
#' 
#' @return A character vector of translated sequences.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(ex.file)
#' translate(fdta$Sequence)
#' 
#' @export
translate <- function( nuc.sequences, M.start=TRUE ){
  nuc.sequences <- gsub( "U", "T", toupper( nuc.sequences ) )
  if( M.start ){
    nuc.sequences <- gsub( "^GTG|^TTG", "ATG", nuc.sequences )
  }
  return( transl( nuc.sequences ) )
}


#' @name reverseComplement
#' @title Reverse-complementation of DNA
#' 
#' @description The standard reverse-complement of nucleotide sequences.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' @param reverse Logical indicating if complement should be reversed.
#' 
#' @details This function uses the Biostrings::reverseComplement function.
#' 
#' @return A character vector of reverse-complemented sequences.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples 
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(ex.file)
#' reverseComplement(fdta$Sequence)
#' 
#' 
#' @export
reverseComplement <- function( nuc.sequences, reverse = TRUE ){
  return(revComp(nuc.sequences, reverse))
}


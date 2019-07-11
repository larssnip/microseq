#' @name translate
#' @title Translation according to the standard genetic code
#' 
#' @description The translation from DNA(RNA) to amino acid sequence according to the standard genetic code.
#' 
#' @param nuc.sequences Character vector containing the nucleotide sequences.
#' @param M.start A logical indicating if the amino acid sequence should start with M regardless of start codon (ATG, GTG or TTG).
#' @param no.stop A logical indicating if terminal stops (*) should be eliminated from the translated sequence
#' @param trans.tab Translation table, either 1 or 4
#' 
#' @details Translating DNA or RNA to protein. Possible start codons are ATG, GTG or TTG. Stop codons are TAA, TGA, TAG. All codons are
#' translated accoring to the Stadard geneti code. This is translation table 1. The only alternative implemented here is translation 
#' table 4, which is used by some bacteria (e.g. Mycoplasma, Mesoplasma). If \code{trans.tab} is 4, the stop codon TGA is 
#' translated to W (Tryptophan).
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
translate <- function(nuc.sequences, M.start = TRUE, no.stop = TRUE, trans.tab = 1){
  nuc.sequences <- gsub("U", "T", toupper(nuc.sequences))
  if(M.start) nuc.sequences <- gsub("^GTG|^TTG", "ATG", nuc.sequences)
  prot.sequences <- transl(nuc.sequences, trans.tab)
  if(no.stop) prot.sequences <- gsub("\\*$", "", prot.sequences)
  return(prot.sequences)
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



#' @name iupac2regex
#' @title Ambiguity symbol conversion
#' @aliases iupac2regex regex2iupac
#' 
#' @description Converting DNA ambiguity symbols to regular expressions, and vice versa.
#' 
#' @usage iupac2regex(sequence)
#' regex2iupac(sequence)
#' 
#' @param sequence Character string containing a DNA sequence.
#' 
#' @details The DNA alphabet may contain ambiguity symbols, e.g. a W means either A or T.
#' When using a regular expression search, these letters must be replaced by the proper
#' regular expression, e.g. W is replaced by [AT] in the string. The \code{iupac2regex} makes this
#' translation, while \code{regex2iupac} cobverts the other way again (replace [AT] with W).
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
#' @export iupac2regex
#' @export regex2iupac
#' 
iupac2regex <- function( sequence ){
  IUPAC <- matrix( c( "W","[AT]",
                      "S","[CG]",
                      "M","[AC]",
                      "K","[GT]",
                      "R","[AG]",
                      "Y","[CT]",
                      "B","[CGT]",
                      "D","[AGT]",
                      "H","[ACT]",
                      "V","[ACG]",
                      "N","[ACGT]" ), ncol=2, byrow=T )
  s <- gsub( "X", "N", toupper( sequence ) )
  for( i in 1:nrow(IUPAC) ){
    s <- gsub( IUPAC[i,1], IUPAC[i,2], s, fixed=T )
  }
  return( s )
}
regex2iupac <- function( sequence ){
  IUPAC <- matrix( c( "W","[AT]",
                      "S","[CG]",
                      "M","[AC]",
                      "K","[GT]",
                      "R","[AG]",
                      "Y","[CT]",
                      "B","[CGT]",
                      "D","[AGT]",
                      "H","[ACT]",
                      "V","[ACG]",
                      "N","[ACGT]" ), ncol=2, byrow=T )
  s <- toupper( sequence )
  for( i in 1:nrow(IUPAC) ){
    s <- gsub( IUPAC[i,2], IUPAC[i,1], s, fixed=T )
  }
  return( s )
}


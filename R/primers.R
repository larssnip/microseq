#' @name amplicon
#' @title Primer matching
#' 
#' @description Extracts subsequences from DNA that matches a specified primer pair.
#' 
#' @param dna Character vector containing the DNA sequences.
#' @param forward String specifying the forward primer.
#' @param reverse String specifying the reverse primer.
#' 
#' @details An amplicon is a subsequence limited by a matching pair of short oligos, called primers.
#' 
#' The \code{forward} primer is a short DNA sequence in the 5' to 3' direction. This can match on both strands of
#' the \code{dna} sequence. The \code{reverse} primer is also a short DNA sequence in 5' to 3' direction, and
#' can also match on both strands of \code{dna}.
#' 
#' For a \code{dna} sequence to produce an amplicon there must be an exact match of the \code{forward}
#' on one strand followed by an exact match of the \code{reverse} on the other strand. The amplicon is
#' the subsequence starting with the \code{forward} and ending with the \code{reverse} primer.
#' 
#' Both primers may contain ambiguity symbols according to the IUPAC standard. 
#' Primers are matched by \code{\link{gregexpr}}, which will not register self-overlapping matches. 
#' In case of multiple (non-overlapping) matches, this function will return all possible amplicons resulting 
#' from the primer matching.
#' 
#' @return A list with the same number of elements as the argument \code{dna}. Each list element
#' contains a string vector with all amplicons resulting from the primer matching. If there is no
#' primer pair match the corresponding string vector is empty.
#' 
#' @author Lars Snipen.
#' 
#' @examples
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(ex.file)
#' amp.lst.1 <- amplicon( fdta$Sequence, forward="AAATTC", reverse="CCAGA" )
#' amp.lst.2 <- amplicon( fdta$Sequence, forward="AANNTC", reverse="CCNGT" ) # more matches due to N's
#' 
#' @export amplicon
#' 
amplicon <- function( dna, forward, reverse ){
  NC <- nchar( dna )
  rvdna <- reverseComplement( dna )
  fwd <- iupac2regex( forward )
  rvs <- iupac2regex( reverse )
  
  ff <- gregexpr( fwd, dna )
  rr <- gregexpr( rvs, rvdna )
  fr <- gregexpr( fwd, rvdna )
  rf <- gregexpr( rvs, dna )
  amplicons <- lapply( 1:length(dna), function(i){
    amps <- character( 0 )
    ss <- ff[[i]]
    ee <- rr[[i]]
    if( ss[1]>0 & ee[1]>0 ){
      ee <- NC[i] - ee + 1
      se <- expand.grid( ss, ee )
      se <- se[which( se[,1]<se[,2] ),]
      if( nrow( se ) >0 ){
        amps <- c( amps,
                   substring( dna[i], se[,1], se[,2] ) )
      }
    }
    ss <- fr[[i]]
    ee <- rf[[i]]
    if( ss[1]>0 & ee[1]>0 ){
      ee <- NC[i] - ee + 1
      se <- expand.grid( ss, ee )
      se <- se[which( se[,1]<se[,2] ),]
      if( nrow( se ) > 0 ){
        amps <- c( amps,
                   reverseComplement( substring( rvdna[i], se[,1], se[,2] ) ) )
      }
    }
    return( amps )
  })
  return( amplicons )
}


#' @name iupac2regex
#' @title Ambiguity symbol conversion
#' @aliases uipac2regex regex2iupac
#' 
#' @description Converting DNA ambiguity symbols to regular expressions, and vice versa.
#' 
#' @usage iupac2regex( s )
#' regex2iupac( s )
#' 
#' @param s Character string containing a DNA sequence.
#' 
#' @details The DNA alphabet may contain ambiguity symbols, e.g. a W means either A or T.
#' When using a regular expression search, these letters must be replaced by the proper
#' regular expression, e.g. W is replaced by [AT] in the string. These function makes this
#' translations
#' 
#' @return A string where the ambiguity symbol has been replaced by a regular expression
#' (\code{iupac2regex}) or a regular expression has been replaced by an ambiguity symbol
#' (\code{regex2iupac}).
#' 
#' @author Lars Snipen.
#' 
#' @examples
#' iupac2regex( "ACWGT" )
#' regex2iupac( "AC[AG]GT" )
#' 
#' @export iupac2regex
#' @export regex2iupac
#' 
iupac2regex <- function( s ){
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
  s <- gsub( "X", "N", toupper( s ) )
  for( i in 1:11 ){
    s <- gsub( IUPAC[i,1], IUPAC[i,2], s, fixed=T )
  }
  return( s )
}
regex2iupac <- function( s ){
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
  s <- toupper( s )
  for( i in 1:11 ){
    s <- gsub( IUPAC[i,2], IUPAC[i,1], s, fixed=T )
  }
  return( s )
}


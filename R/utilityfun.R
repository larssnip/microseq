#' @name gregexpr
#' @title Extended \code{\link{gregexpr}} with substring retrieval
#' 
#' @description An extension of the function \code{base::gregexpr} enabling retrieval of the 
#' matching substrings.
#' 
#' @param pattern Character string containing a \link{regular expression} (or character string for 
#' \code{fixed = TRUE}) to be matched in the given character vector.  Coerced by \code{\link{as.character}} 
#' to a character string if possible.  If a character vector of length 2 or more is supplied, the first element
#' is used with a warning.  Missing values are not allowed.
#' @param text A character vector where matches are sought, or an object which can be coerced by
#' \code{\link{as.character}} to a character vector.
#' @param ignore.case If \code{FALSE}, the pattern matching is \emph{case sensitive} and if \code{TRUE},
#' case is ignored during matching.
#' @param perl Logical. Should perl-compatible regexps be used? Has priority over \code{extended}.
#' @param fixed Logical. If \code{TRUE}, \samp{pattern} is a string to be matched as is. Overrides 
#' all conflicting arguments.
#' @param useBytes Logical. If \code{TRUE} the matching is done byte-by-byte rather than character-by-character.
#' See \code{base::gregexpr} for details.
#' @param extract Logical indicating if matching substrings should be extracted and returned.
#' 
#' @details Extended version of \code{base:gregexpr} that enables the return of the substrings matching
#' the pattern. The last argument \samp{extract} is the only difference to \code{base::gregexpr}. The default
#' behaviour is identical to \code{base::gregexpr}, but setting \code{extract=TRUE} means the matching substrings
#' are returned.
#' 
#' @return It will either return what the \code{base::gregexpr} would (\code{extract=FALSE}) or a \samp{list}
#' of substrings matching the pattern (\code{extract=TRUE}). There is one \samp{list} element for each string in
#' \samp{text}, and each list element contains a character vector of all matching substrings in the corresponding
#' entry of \samp{text}.
#' 
#' @author Lars Snipen and Kristian Liland.
#' 
#' @seealso \code{\link[base]{gregexpr}}
#' 
#' @examples
#' sequences <- c("ACATGTCATGTCC", "CTTGTATGCTG")
#' gregexpr("ATG", sequences, extract = TRUE)
#' 
#' @keywords gregexpr
#' 
#' @export gregexpr
#' 
gregexpr <- function(pattern, text, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, extract = FALSE){
  lst <- base::gregexpr(pattern, text, ignore.case, perl, fixed, useBytes)
  if(extract){
    lst <- lapply(1:length(lst), function(i){substring(text[i], lst[[i]], lst[[i]]+attr(lst[[i]], "match.length") - 1)})
    }
  return(lst)
}
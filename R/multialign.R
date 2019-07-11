#' @name msalign
#' @title Multiple alignment
#' 
#' @description Quickly computing a smallish multiple sequence alignment.
#' 
#' @param fdta A \code{Fasta} object with input sequences.
#' @param machine Function that does the 'dirty work'.
#' 
#' @details This function computes a multiple sequence alignment given a set of sequences in a \code{Fasta} object,
#' see \code{\link{readFasta}} for more on \code{Fasta} objects.
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
#' @return Results are returned as a \code{Fasta} object.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{muscle}}, \code{\link{msaTrim}}.
#' 
#' @importFrom stats rnorm
#' 
#' @examples
#' \dontrun{
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(ex.file)
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



#' @name muscle
#' @title Multiple alignment using MUSCLE
#' 
#' @description Computing a multiple sequence alignment using the MUSCLE software.
#' 
#' @param in.file Name of FASTA-file with input sequences.
#' @param out.file Name of file to store the result.
#' @param quiet Logical, \code{quiet=FALSE} produces screen output during computations.
#' @param diags Logical, \code{diags=TRUE} gives faster but less reliable alignment.
#' @param maxiters Maximum number of iterations.
#' 
#' @details The software MUSCLE (Edgar, 2004) must be installed and available on the system. Test this by typing
#' \code{system("muscle")} in the Console, and some sensible output should be produced. NOTE: The executable must
#' be named \code{muscle} on your system, no version numbers etc. For more details
#' on MUSCLE, see \url{http://www.drive5.com/muscle}.
#' 
#' By default \code{diags=FALSE} but can be set to \code{TRUE} to increase speed. This should be done only 
#' if sequences are highly similar.
#' 
#' By default \code{maxiters=16}. If you have a large number of sequences (a few thousand), or they are 
#' very long, then this may be too slow for practical use. A good compromise between speed and accuracy is to run
#' just the first two iterations of the algorithm. On average, this gives accuracy equal to T-Coffee
#' and speeds much faster than CLUSTALW. This is done by the option \code{maxiters=2}.
#' 
#' @return The result is written to the file specified in \code{out.file}.
#' 
#' @author Lars Snipen.
#' 
#' @references Edgar, R.C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high 
#' throughput, Nucleic Acids Res, 32, 1792-1797.
#' 
#' @seealso \code{\link{msaTrim}}.
#' 
#' @examples
#' \dontrun{
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' muscle(in.file=ex.file,out.file="deleteMe.fasta")
#' }
#' 
#' @export muscle
#' 
muscle <- function(in.file, out.file, quiet = FALSE, diags = FALSE, maxiters = 16){
  if(available.external("muscle")){
    dtxt <- ifelse(diags, "-diags", "")
    itxt <- paste("-maxiters", maxiters)
    qt <- ifelse(quiet, "-quiet", "")
    cmd <- paste("muscle", dtxt, itxt, qt, "-in", in.file, "-out", out.file)
    system(cmd)
  }
}



#' @name cmalign
#' @title Multiple alignment using Infernal
#' 
#' @description Computing a multiple sequence alignment using the Infernal software.
#' 
#' @param in.file Name of FASTA-file with input sequences.
#' @param out.file Name of file to store the result.
#' @param CM.file Name of file with correlation model.
#' @param threads Number of CPU's to use
#' 
#' @details The software Infernal (Nawrocki&Eddy, 2013) must be installed and available on the system. Test 
#' this by typing \code{system("cmalign -h")} in the Console, and some sensible output should be produced. 
#' For more details on Infernal, see http://eddylab.org/infernal/.
#' 
#' This function is most typically used to align 16S rRNA sequences.
#' 
#' The \code{cmalign} function will produce a multiple alignment, like e.g. \code{muscle}, but makes use
#' of a \emph{correlation model} to do so. A correlation model means in this case a description of
#' how various bases have a long-range relation, due to folding of the sequence. This means that you can
#' only use this function to align sequences for which you have such correlation models. Such models
#' are typically available for a number of RNA-families, see below. 
#' 
#' The argument \code{CM.file} is the name of a file with a valid correlation model, e.g. one downloaded
#' from the \href{http://rfam.xfam.org/}{Rfam database}. See examples below for the 16S model supplied with this
#' package.
#' 
#' @return The result is written to the file specified in \code{out.file}.
#' 
#' @author Lars Snipen.
#' 
#' @references E.P. Nawrocki and S.R. Eddy,  Infernal 1.1: 100-fold faster RNA homology searches, 
#' Bioinformatics 29:2933-2935 (2013). 
#' 
#' @seealso \code{\link{msaTrim}}.
#' 
#' @examples
#' \dontrun{
#' in.file <- file.path(file.path(path.package("microseq"),"extdata"),"16S.fasta")
#' cm.file <- file.path(file.path(path.package("microseq"),"extdata"),"ssu_bacteria.cm")
#' cmalign(in.file,"msa_infernal.fasta",cm.file)
#' }
#' 
#' 
#' @export cmalign
#' 
cmalign <- function(in.file, out.file, CM.file, threads = 1){
  if(available.external('infernal')){
    cmd <- paste("cmalign --cpu ",
                 threads,
                 "-o cmalign.stk",
                 "--outformat Pfam",
                 CM.file,
                 in.file)
    system(cmd)
    lines <- readLines("cmalign.stk")
    fdta <- stk2fasta(lines)
    writeFasta(fdta, out.file)
    file.remove("cmalign.stk")
  }
}

stk2fasta <- function(lines){
  lines <- lines[nchar(lines)>0]
  lines <- lines[lines != "//"]
  lines <- lines[grep( "^#", lines, invert = T )]
  lst <- strsplit(lines, split = " ")
  header <- sapply(lst, function(x){x[1]})
  sequence <- toupper(sapply(lst, function(x){x[length(x)]}))
  
  fdta <- data.frame(Header = header,
                     Sequence = gsub("U", "T", gsub("\\.", "-", sequence)),
                     stringsAsFactors = F)
  class(fdta) <- c("Fasta", "data.frame")
  return(fdta)
}



#' @name msaTrim
#' @title Trimming multiple sequence alignments
#' 
#' @description Trimming a multiple sequence alignment by discarding columns with too many gaps.
#' 
#' @param msa A \code{Fasta} object containing a multiple alignment.
#' @param gap.end Fraction of gaps tolerated at the ends of the alignment (0-1).
#' @param gap.mid Fraction of gaps tolerated inside the alignment (0-1).
#' 
#' @details A multiple alignment is trimmed by removing columns with too many indels (gap-symbols). Any 
#' columns containing a fraction of gaps larger than \code{gap.mid} are discarded. For this reason, \code{gap.mid}
#' should always be farily close to 1.0 therwise too many columns may be discarded, destroying the alignment.
#' 
#' Due to the
#' heuristics of multiple alignment methods, both ends of the alignment tend to be uncertain and most
#' of the trimming should be done at the ends. Starting at each end, columns are discarded as long as their fraction of gaps
#' surpasses \code{gap.end}. Typically \code{gap.end} can be much smaller than \code{gap.mid}, but if 
#' set too low you risk that all columns are discarded!
#' 
#' @return The trimmed alignment is returned as a \code{Fasta} object.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{muscle}}, \code{\link{cmalign}}.
#' 
#' @examples 
#' \dontrun{
#' msa.file <- file.path(file.path(path.package("microseq"),"extdata"),"msa.fasta")
#' msa <- readFasta(msa.file)
#' print(nchar(msa$Sequence))
#' msa.trimmed <- msaTrim(msa)
#' print(nchar(msa.trimmed$Sequence))
#' }
#' 
#' @export msaTrim
#' 
msaTrim <- function(msa, gap.end = 0.5, gap.mid = 0.9){
  nc <- unique(nchar(msa$Sequence))
  if(length(nc) != 1) stop("This is not a multiple alignment, sequences have different lengths!")
  cmat <- t(sapply(strsplit(msa$Sequence, split = ""), function(x){x}))
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


## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(what){
  if(what == "muscle"){
    chr <- NULL
    try(chr <- system('muscle -version', intern = TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('MUSCLE was not found by R.',
                 'Please install MUSCLE from: http://www.drive5.com/muscle/downloads.htm.',
                 'After installation, make sure MUSCLE can be run from a shell/terminal, ',
                 'using the command \'muscle\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  
  if(what == "infernal"){
    chr <- NULL
    try(chr <- system('cmalign -h', intern = TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('Infernal was not found by R.',
                 'Please install Infernal from: http://eddylab.org/infernal/.',
                 'After installation, make sure Infernal can be run from a shell/terminal, ',
                 'using the command \'cmalign\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

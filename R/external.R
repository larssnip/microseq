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


#' @name findrRNA
#' @title Finding rRNA genes
#' 
#' @description Locating all rRNA genes in genomic DNA using the barrnap software.
#' 
#' @param genome.file A fasta-formatted file with the genome sequence(s).
#' @param bacteria Logical, the genome is either a bacteria (default) or an archea.
#' @param cpu Number of CPUs to use, default is 1.
#' 
#' @details The external software barrnap is used to scan through a prokaryotic genome to detect the
#' rRNA genes (5S, 16S, 23S). This free software can be installed from https://github.com/tseemann/barrnap.
#' 
#' @return A \code{gff.table} (see \code{\link{readGFF}} for details) with one row for each detected
#' rRNA sequence.
#' 
#' @note The barrnap software must be installed on the system for this function to work, i.e. the command
#' \samp{system("barrnap --help")} must be recognized as a valid command if you run it in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}.
#' 
#' @examples
#' \dontrun{
#' # This example requires the external barrnap software
#' # Using a genome file in this package.
#' xpth <- file.path(path.package("microseq"),"extdata")
#' genome.file <- file.path(xpth,"small_genome.fasta")
#' 
#' # Searching for rRNA sequences, and inspecting
#' gff.tbl <- findrRNA(genome.file)
#' print(gff.table)
#' 
#' # Retrieving the sequences
#' genome <- readFasta(genome.file)
#' rRNA <- gff2fasta(gff.tbl, genome)
#' }
#' 
#' @export findrRNA
#' 
findrRNA <- function(genome.file, bacteria = TRUE, cpu = 1){
  if(available.external("barrnap")){
    kingdom <- ifelse(bacteria, "bac", "arc")
    tmp.file <- tempfile(pattern = "barrnap", fileext = ".gff")
    system(paste("barrnap --quiet --kingdom", kingdom, genome.file, ">", tmp.file))
    gff.table <- readGFF(tmp.file)
    file.remove(tmp.file)
    return(gff.table)
  }
}










## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(what){
  if(what == "muscle"){
    chr <- NULL
    try(chr <- system2("muscle", args = "-version", stdout = TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('muscle was not found by R.',
                 'Please install muscle from: http://www.drive5.com/muscle/downloads.htm.',
                 'After installation, make sure muscle can be run from a shell/terminal, ',
                 'using the command \'muscle -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  
  if(what == "barrnap"){
    chr <- NULL
    try(chr <- system2("barrnap", args = "--version", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('barrnap was not found by R.',
                 'Please install barrnap from: https://github.com/tseemann/barrnap',
                 'After installation, make sure barrnap can be run from a shell/terminal, ',
                 'using the command \'barrnap --help\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

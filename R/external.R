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
#' @description Finding rRNA genes in genomic DNA using the barrnap software.
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



#' @name findGenes
#' @title Finding coding genes
#' 
#' @description Finding coding genes in genomic DNA using the Prodigal software.
#' 
#' @param genome.file A fasta-formatted file with the genome sequence(s).
#' @param faa.file If provided, prodigal will output all proteins to this fasta-file (text).
#' @param ffn.file If provided, prodigal will output all DNA sequences to this fasta-file (text).
#' @param proc Either \code{"single"} or \code{"meta"}, see below.
#' @param trans.tab Either 11 or 4 (see below).
#' @param mask.N Turn on masking of N's (logical)
#' @param bypass.SD Bypass Shine-Dalgarno filter (logical)
#' 
#' @details The external software Prodigal is used to scan through a prokaryotic genome to detect the protein
#' coding genes. This free software can be installed from https://github.com/hyattpd/Prodigal.
#' 
#' In addition to the standard output from this function, fasta-files with protein and/or DNA sequences may
#' be produced directly by providing filenames in \code{faa.file} and \code{ffn.file}.
#' 
#' The input \code{proc} allows you to såecify i the input data should be treated as a single genome
#' (default) or as a metagenome.
#' 
#' The translation table is by default 11 (the standard code), but table 4 should be used for Mycoplasma etc.
#' 
#' The \code{mask.N} will prevent genes having runs of N inside. The \code{bypass.SD} turn off the search
#' for a Shine-Dalgarno motif.
#' 
#' @return A GFF-table (see \code{\link{readGFF}} for details) with one row for each detected
#' coding gene.
#' 
#' @note The prodigal software must be installed on the system for this function to work, i.e. the command
#' \samp{system("prodigal -h")} must be recognized as a valid command if you run it in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}.
#' 
#' @examples
#' \dontrun{
#' # This example requires the external prodigal software
#' # Using a genome file in this package.
#' xpth <- file.path(path.package("microseq"),"extdata")
#' genome.file <- file.path(xpth,"small_genome.fasta")
#' 
#' # Searching for coding sequences, and inspecting
#' gff.tbl <- findGenes(genome.file)
#' 
#' # Retrieving the sequences
#' genome <- readFasta(genome.file)
#' cds <- gff2fasta(gff.tbl, genome)
#' }
#' 
#' @export findGenes
#' 
findGenes <- function(genome.file, faa.file = "", ffn.file = "", proc = "single",
                      trans.tab = 11, mask.N = FALSE, bypass.SD = FALSE){
  if(available.external("prodigal")){
    if(nchar(faa.file) > 0) faa.file <- paste("-a", faa.file)
    if(nchar(ffn.file) > 0) ffn.file <- paste("-d", ffn.file)
    mask.N <- ifelse(mask.N, "-m", "")
    bypass.SD <- ifelse(bypass.SD, "-n", "")
    tmp.file <- tempfile(pattern = "prodigal", fileext = "gff")
    command <- paste("prodigal -q -f gff",
                     protein.file,
                     mrna.file,
                     "-p", proc,
                     mask.N,
                     bypass.SD,
                     "-g", trans.tab,
                     "-i", genome.file,
                     "-o", tmp.file)
    system(command)
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
  
  if(what == "prodigal"){
    chr <- NULL
    try(chr <- system2("prodigal", args = "-v", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('prodigal was not found by R.',
                 'Please install barrnap from: https://github.com/hyattpd/Prodigal',
                 'After installation, make sure prodigal can be run from a shell/terminal, ',
                 'using the command \'prodigal -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}
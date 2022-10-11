#' @name muscle
#' @title Multiple alignment using MUSCLE
#' 
#' @description Computing a multiple sequence alignment using the MUSCLE software.
#' 
#' @param in.file Name of FASTA file with input sequences.
#' @param out.file Name of file to store the result.
#' @param muscle.exe Command to run the external software muscle on the system (text).
#' @param quiet Logical, \code{quiet = FALSE} produces screen output during computations.
#' @param diags Logical, \code{diags = TRUE} gives faster but less reliable alignment.
#' @param maxiters Maximum number of iterations.
#' 
#' @details The software MUSCLE (Edgar, 2004) must be installed and available on the system. The text in
#' \code{muscle.exe} must contain the exact command to invoke muscle on the system.
#' 
#' By default \code{diags = FALSE} but can be set to \code{TRUE} to increase speed. This should be done
#' only if sequences are highly similar.
#' 
#' By default \code{maxiters = 16}. If you have a large number of sequences (a few thousand), or they are 
#' very long, then this may be too slow for practical use. A good compromise between speed and accuracy
#' is to run just the first two iterations of the algorithm. On average, this gives accuracy equal to
#' T-Coffee and speeds much faster than CLUSTALW. This is done by the option \code{maxiters = 2}.
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
#' fa.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.faa")
#' muscle(in.file = fa.file, out.file = "delete_me.msa")
#' }
#' 
#' @export muscle
#' 
muscle <- function(in.file, out.file, muscle.exe = "muscle", quiet = FALSE, diags = FALSE, maxiters = 16){
  if(available.external(muscle.exe)){
    dtxt <- ifelse(diags, "-diags", "")
    itxt <- paste("-maxiters", maxiters)
    qt <- ifelse(quiet, "-quiet", "")
    cmd <- paste(mucle.exe, dtxt, itxt, qt, "-in", in.file, "-out", out.file)
    system(cmd)
  }
}


#' @name findrRNA
#' @title Finding rRNA genes
#' 
#' @description Finding rRNA genes in genomic DNA using the barrnap software.
#' 
#' @param genome A table with columns Header and Sequence, containing the genome sequence(s).
#' @param barrnap.exe Command to run the external software barrnap on the system (text).
#' @param bacteria Logical, the genome is either a bacteria (default) or an archea.
#' @param cpu Number of CPUs to use, default is 1.
#' 
#' @details The external software barrnap is used to scan through a prokaryotic genome to detect the
#' rRNA genes (5S, 16S, 23S).
#' The text in \code{barrnap.exe} must contain the exact command to invoke barrnap on the system.
#' 
#' @return A GFF-table (see \code{\link{readGFF}} for details) with one row for each detected
#' rRNA sequence.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}.
#' 
#' @examples
#' \dontrun{
#' # This example requires the external barrnap software
#' # Using a genome file in this package.
#' genome.file <- file.path(path.package("microseq"),"extdata","small.fna")
#' 
#' # Searching for rRNA sequences, and inspecting
#' genome <- readFasta(genome.file)
#' gff.tbl <- findrRNA(genome)
#' print(gff.table)
#' 
#' # Retrieving the sequences
#' rRNA <- gff2fasta(gff.tbl, genome)
#' }
#' 
#' @export findrRNA
#' 
findrRNA <- function(genome, barrnap.exe = "barrnap", bacteria = TRUE, cpu = 1){
  if(available.external(barrnap.exe)){
    if(sum(c("Header", "Sequence") %in% colnames(genome)) != 2) stop("First argument must be table with columns Header and Sequence")
    if(nrow(genome) == 0) stop("Genome is empty (no sequences)")
    kingdom <- ifelse(bacteria, "bac", "arc")
    tmp.fna <- tempfile(pattern = "genome", fileext = ".fna")
    writeFasta(genome, out.file = tmp.fna)
    tmp.gff <- tempfile(pattern = "barrnap", fileext = ".gff")
    system(paste("barrnap --quiet --kingdom", kingdom, tmp.fna, ">", tmp.gff))
    gff.table <- readGFF(tmp.gff)
    ok <- file.remove(tmp.fna, tmp.gff)
    return(gff.table)
  }
}



#' @name findGenes
#' @title Finding coding genes
#' 
#' @description Finding coding genes in genomic DNA using the Prodigal software.
#' 
#' @param genome A table with columns Header and Sequence, containing the genome sequence(s).
#' @param faa.file If provided, prodigal will output all proteins to this fasta-file (text).
#' @param ffn.file If provided, prodigal will output all DNA sequences to this fasta-file (text).
#' @param prodigal.exe Command to run the external software prodigal on the system (text).
#' @param proc Either \code{"single"} or \code{"meta"}, see below.
#' @param trans.tab Either 11 or 4 (see below).
#' @param mask.N Turn on masking of N's (logical)
#' @param bypass.SD Bypass Shine-Dalgarno filter (logical)
#' 
#' @details The external software Prodigal is used to scan through a prokaryotic genome to detect the protein
#' coding genes. The text in \code{prodigal.exe} must contain the exact command to invoke barrnap on the system.
#' 
#' In addition to the standard output from this function, FASTA files with protein and/or DNA sequences may
#' be produced directly by providing filenames in \code{faa.file} and \code{ffn.file}.
#' 
#' The input \code{proc} allows you to specify if the input data should be treated as a single genome
#' (default) or as a metagenome. In the latter case the \code{genome} are (un-binned) contigs.
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
#' genome.file <- file.path(path.package("microseq"),"extdata","small.fna")
#' 
#' # Searching for coding sequences, this is Mycoplasma (trans.tab = 4)
#' genome <- readFasta(genome.file)
#' gff.tbl <- findGenes(genome, trans.tab = 4)
#' 
#' # Retrieving the sequences
#' cds.tbl <- gff2fasta(gff.tbl, genome)
#' 
#' # You may use the pipe operator
#' library(ggplot2)
#' readFasta(genome.file) %>% 
#'   findGenes(trans.tab = 4) %>% 
#'   filter(Score >= 50) %>% 
#'   ggplot() +
#'   geom_histogram(aes(x = Score), bins = 25)
#' }
#' 
#' @export findGenes
#' 
findGenes <- function(genome, prodigal.exe = "prodigal", faa.file = "", ffn.file = "", proc = "single",
                      trans.tab = 11, mask.N = FALSE, bypass.SD = FALSE){
  if(available.external(prodigal.exe)){
    if(sum(c("Header", "Sequence") %in% colnames(genome)) != 2) stop("First argument must be table with columns Header and Sequence")
    if(nrow(genome) == 0) stop("Genome is empty (no sequences)")
    if(nchar(faa.file) > 0) faa.file <- paste("-a", faa.file)
    if(nchar(ffn.file) > 0) ffn.file <- paste("-d", ffn.file)
    mask.N <- ifelse(mask.N, "-m", "")
    bypass.SD <- ifelse(bypass.SD, "-n", "")
    tmp.fna <- tempfile(pattern = "genome", fileext = "fna")
    writeFasta(genome, out.file = tmp.fna)
    tmp.gff <- tempfile(pattern = "prodigal", fileext = "gff")
    command <- paste("prodigal -q -f gff",
                     faa.file,
                     ffn.file,
                     "-p", proc,
                     mask.N,
                     bypass.SD,
                     "-g", trans.tab,
                     "-i", tmp.fna,
                     "-o", tmp.gff)
    system(command)
    gff.table <- readGFF(tmp.gff)
    file.remove(tmp.fna, tmp.gff)
    return(gff.table)
  }
}


## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(exe){
  chr <- NULL
  ok <- try(chr <- system2(exe, stdout = TRUE), silent = TRUE)
  if(length(grep("Error", ok[1])) > 0){
    stop("Cannot run ", exe, " from R")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

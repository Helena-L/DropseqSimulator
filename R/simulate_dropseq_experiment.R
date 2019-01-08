#' Simulator: simulate_dropseq_experiment
#' Simulate Dropseq RNA-seq reads
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, given gene number and cell number
#' @param fasta path to FASTA file containing transcripts from which to simulate reads.
#' @param ngenes number of genes to simulate
#' @param ncells number of cells to simulate
#' @param polyAnum minimum number of 'A's in a polyA region, integer. A region should contain at least n continous 'A's to be considered as a polyA region.
#' @param outdir path to folder where simulated reads should be written. By default, reads written to the working directory.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export

simulate_dropseq_experiment = function(fasta=NULL, ngenes=NULL, ncells=NULL, polyAnum=15, outdir='.') {
  if(!is.null(fasta)){
    transcripts = readDNAStringSet(fasta)
  }
  else{
    stop('must provide fasta file')
  }

  if(is.null(ngenes) || is.null(ncells)){
    stop('must provide # genes and # cells')
  }

  stopifnot(ngenes > 0)
  stopifnot(ncells > 0)

  stopifnot(ngenes == length(transcripts))

  readmat = generate_countmat(ngenes, ncells)

  simulate_dropseq_experiment_countmat(fasta=fasta, readmat=readmat, polyAnum=polyAnum, outdir=outdir)
}

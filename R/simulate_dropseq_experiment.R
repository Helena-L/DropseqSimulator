#' Simulator: simulate_dropseq_experiment
#' Simulate Dropseq RNA-seq reads
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, given gene number and cell number
#' @param fasta path to FASTA file containing transcripts from which to simulate reads.
#' @param ngenes number of genes to simulate
#' @param ncells number of cells to simulate
#' @param libloc location parameter for the library size log-normal distribution, together with \code{libscale} controls library size. Determine different sequencing depth (number of reads) between cells.
#' @param libscale scale parameter for the library size log-normal distribution, together with \code{libloc} controls library size. Determine different sequencing depth (number of reads) between ceels.
#' @param polyAnum minimum number of 'A's in a polyA region, integer. A region should contain at least n continous 'A's to be considered as a polyA region.
#' @param bias polyA sampling bias, one of 'empirical', 'custom' (default 'empirical'). If 'empirical', the built-in model is adopted.
#' This model is trained on data sets from the Drop-seq experiments described in the Cell paper (Macosko et al, 2015).
#' If 'custom', polyA sampling bias is captured from user input model. The path of user input model is specified in \code{model}.
#' @param model path to polyA sampling bias model if \code{bias} set to 'custom'. Model format should be a data frame with 2 columns.
#' Column 1 is distance to polyA region, while Column 2 is probability of getting a fragment at that position. Column 2 sums to 1.
#' @param outdir path to folder where simulated reads should be written. By default, reads written to the working directory.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export

simulate_dropseq_experiment = function(fasta=NULL, ngenes=NULL, ncells=NULL, libloc=5, libscale=0.2, polyAnum=15, bias='empirical', model=NULL, outdir='.') {
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

  readmat = generate_countmat(ngenes, ncells, libloc, libscale)

  simulate_dropseq_experiment_countmat(fasta=fasta, readmat=readmat, polyAnum=polyAnum, outdir=outdir, bias=bias, model=model)
}

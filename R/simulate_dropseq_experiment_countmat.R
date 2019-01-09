#' Simulator: simulate_dropseq_experiment_countmat
#' Simulate Dropseq RNA-seq reads
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, given a countmat matrix
#' @param fasta path to FASTA file containing transcripts from which to simulate reads.
#' @param readmat a gene-cell matrix, each entry represents number of reads to simulate.
#' @param polyAnum minimum number of 'A's in a polyA region, integer. A region should contain at least n continous 'A's to be considered as a polyA region.
#' @param bias polyA sampling bias, one of 'empirical', 'custom', 'naive' (default 'empirical'). If 'empirical', the built-in model is adopted.
#' This model is trained on data sets from the Drop-seq experiments described in the Cell paper (Macosko et al, 2015).
#' If 'custom', polyA sampling bias is captured from user input model. The path of user input model is specified in \code{model}.
#' If 'naive', sampling probability is calculated as (1/distance to polyA region).
#' @param model path to polyA sampling bias model if \code{bias} set to 'custom'. Model format should be a data frame with 2 columns.
#' Column 1 is distance to polyA region, while Column 2 is probability of getting a fragment at that position. Column 2 sums to 1.
#' @param outdir path to folder where simulated reads should be written. By default, reads written to the working directory.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export

simulate_dropseq_experiment_countmat = function(fasta=NULL, readmat, polyAnum=15, bias='empirical', model=NULL, outdir='.') {
  if(!is.null(fasta)){
    transcripts = readDNAStringSet(fasta)
  }
  else{
    stop('must provide fasta file')
  }

  stopifnot(class(readmat) == 'matrix')
  stopifnot(nrow(readmat) == length(transcripts))

  import_path = system.file("polyester", "R", package="Simulator")
  import_files = list.files(import_path)
  for(i in import_files){
    source(paste(import_path, '/', i, sep=''))
  }

  if(!(bias %in% c('empirical', 'custom', 'naive'))){
    stop('bias should be empirical or custom or naive')
  }
  if(bias == 'custom'){
    if(is.null(model)){
      stop('need to provide polyA bias model')
    }
  }
  # simulate sequence reads
  message("Simulating sequence reads...")
  simulate_experiment_countmat(fasta, readmat=readmat, outdir=outdir, paired=FALSE, readlen=50, fraglen=100, fragsd=10, bias='dropseqf_polyA', polyAnum=polyAnum, polyAbias=bias, polyAmodel=model)
  # simulate barcodes
  message("Simulating barcodes...")
  generate_barcodes(outdir)
  message("Done!")
}

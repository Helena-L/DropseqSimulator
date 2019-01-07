#' Simulator: simulate_dropseq_experiment_countmat
#' Simulate Dropseq RNA-seq reads
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, given a countmat matrix
#' @param fasta path to FASTA file containing transcripts from which to simulate reads.
#' @param readmat a gene-cell matrix, each entry represents number of reads to simulate
#' @param outdir path to folder where simulated reads should be written. By default, reads written to the working directory.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export

simulate_dropseq_experiment_countmat = function(fasta=NULL, readmat, outdir='.') {
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
  simulate_experiment_countmat(fasta, readmat=readmat, outdir=outdir, paired=FALSE, readlen=50, fraglen=100, fragsd=10, bias='dropseqf_polyA')
}

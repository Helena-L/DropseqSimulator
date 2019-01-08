#' generate barcodes
#' generate cell barcode and molecolar barcode for each read
#' @param dirpath file path specifying where sequencing reads should be written
#' @return No return, but simulated reads are written to \code{filepath}.
#' @export

generate_barcodes = function(dirpath) {
  files = list.files(dirpath)
  ncells = length(files)

  cell_barcode_pool = c()
  i = 1
  while(i <= ncells){
    barcode = ''
    for(t in 1:12){
      base = sample(c('A', 'T', 'C', 'G'),1)
      barcode = paste(barcode, base, sep="")
    }
    if(!is.element(barcode, cell_barcode_pool)){
      cell_barcode_pool = c(cell_barcode_pool, barcode)
      i = i+1
    }
  }

  for(i in 1:ncells){
    if(!grepl("fasta", files[i])){
      next
    }
    filepath = paste(dirpath, "/", files[i], sep="")
    fasta = readDNAStringSet(filepath)
    nreads = length(fasta)
    molecular_barcode_pool = c()
    j = 1
    while(j <= nreads){
      barcode = ''
      for(t in 1:8){
        base = sample(c('A', 'T', 'C', 'G'),1)
        barcode = paste(barcode, base, sep="")
      }
      if(!is.element(barcode, molecular_barcode_pool)){
        molecular_barcode_pool = c(molecular_barcode_pool, barcode)
        j = j+1
      }
    }
    # add cell barcode and molecular barcode together
    cm_barcodes = paste(cell_barcode_pool[i], molecular_barcode_pool, sep="")
    # write cm_barcodes to fasta
    filename = paste(substr(files[i], 1, width(files[i])-6), '_barcodes.fasta', sep="")
    cm_barcodes_fasta = DNAStringSet(cm_barcodes)
    names(cm_barcodes_fasta) = names(fasta)
    writeXStringSet(cm_barcodes_fasta, paste(dirpath, "/", filename, sep=""), format="fasta")
  }
}

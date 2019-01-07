#' find reverse complement sequence given a transcript
#'
#' @param transcript transcript sequence, DNAString.
#' @return reverse complement sequence of input transcript, DNAString.
#' @export

find_reverse_complement = function(transcript) {
  seq = as.character(transcript)
  seq = strsplit(seq, "")[[1]]
  rc_seq = c()
  for(ch in seq){
    if(ch == 'A'){
      add_ch = 'T'
    }
    else if(ch == 'T'){
      add_ch = 'A'
    }
    else if(ch == 'C'){
      add_ch = 'G'
    }
    else{
      add_ch = 'C'
    }
    rc_seq = c(rc_seq, add_ch)
  }
  rc_seq = rev(rc_seq)
  rc_seq = paste(rc_seq,collapse="")
  rc_seq = DNAString(rc_seq)
  return (rc_seq)
}

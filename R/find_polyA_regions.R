#' find polyA regions given a transcript
#'
#' @param transcript transcript sequence, DNAString.
#' @return polyA regions along the transcript, list.
#' @export

find_polyA_regions = function(transcript) {
  polyAs = "AAAAAAAAAAAAAAA"
  seq = as.character(transcript)
  polyA_regions = list()
  num = 0
  if(grepl(polyAs, seq)){
    polyA_st_pos = regexpr(polyAs, seq)[1]
    polyA_end_pos = polyA_st_pos
    while(polyA_end_pos < width(seq)){
      while(polyA_end_pos <= width(seq) && substr(seq, polyA_end_pos, polyA_end_pos) == 'A'){
        polyA_end_pos = polyA_end_pos+1
      }
      num = num+1
      polyA_regions[[num]] = c(polyA_st_pos, polyA_end_pos-1)
      polyA_st_pos_sub = regexpr(polyAs, substr(seq, polyA_end_pos, width(seq)))[1]
      if(polyA_st_pos_sub == -1){
        break
      }
      polyA_st_pos = polyA_st_pos_sub+polyA_end_pos-1
      polyA_end_pos = polyA_st_pos
    }
  }
  return (polyA_regions)
}

#' given a transcript and fragment length, calculate probability of generating reads along the transcript
#'
#' @param transcript transcript sequence, DNAString.
#' @param fraglen fragment length, integer.
#' @return probability of generating reads along the transcript
#' @export
calculate_sample_prob_naive = function(transcript_len, fraglen) {
  prob_list = c()
  end_pos = transcript_len-fraglen+1
  for(i in 1:end_pos){
    prob_list = c(prob_list, 1/i)
  }
  prob_list = rev(prob_list)
  prob_list_norm = c()
  for(prob in prob_list){
    prob_list_norm = c(prob_list_norm, prob / sum(prob_list))
  }
  return (prob_list_norm)
}

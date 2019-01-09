#' given a transcript and fragment length, calculate probability of generating reads along the transcript
#'
#' @param transcript transcript sequence, DNAString.
#' @param fraglen fragment length, integer.
#' @param polyAnum minimum number of 'A's in a polyA region, integer. A region should contain at least n continous 'A's to be considered as a polyA region.
#' @param bias polyA sampling bias, one of 'empirical', 'custom', 'naive' (default 'empirical'). If 'empirical', the built-in model is adopted.
#' This model is trained on data sets from the Drop-seq experiments described in the Cell paper (Macosko et al, 2015).
#' If 'custom', polyA sampling bias is captured from user input model. The path of user input model is specified in \code{model}.
#' If 'naive', sampling probability is calculated as (1/distance to polyA region).
#' @param model path to polyA sampling bias model if \code{bias} set to 'custom'. Model format should be a data frame with 2 columns.
#' Column 1 is distance to polyA region, while Column 2 is probability of getting a fragment at that position. Column 2 sums to 1.
#' @return probability of generating reads along the transcript
#' @export
calculate_sample_prob_polyA = function(transcript, fraglen, polyAnum=15, bias='empirical', model=NULL) {
  rc_transcript = find_reverse_complement(transcript)
  # find the polyA region in the input transcript
  polyA_regions = c(find_polyA_regions(transcript, polyAnum), find_polyA_regions(rc_transcript, polyAnum))
  # for each base, calculate the distance to nearest polyA region
  end_pos = length(transcript)-fraglen+1
  transcript_polyA_arr = c()
  for(i in 1:end_pos){
    # calculate the distance of i to its nearest polyA region
    x_distance_to_polyA = length(transcript)
    for(region in polyA_regions){
      start = region[1]
      end = region[2]
      if(i <= end){
        if((start <= i) && (i <= end)){
          x_distance_to_polyA = 0
          break
        }
        x_distance_to_polyA = min(x_distance_to_polyA, start-i)
      }
    }
    # compare the distance to polyA to distance to 3' end
    x_distance_to_polyA = min(x_distance_to_polyA, length(transcript)-i)
    transcript_polyA_arr = c(transcript_polyA_arr, x_distance_to_polyA)
  }

  if(bias == 'empirical'){
    # read in trained model:  distance to polyA, num of reads
    file_path = system.file("extdata", "polyA_model.Rda", package="Simulator")
    load(file_path)
    read_num_list = c()
    for(x in transcript_polyA_arr){
      read_num_list = c(read_num_list, polyA_model$number[x+1])
    }
    prob_list_norm = c()
    for(read_num in read_num_list){
      prob_list_norm = c(prob_list_norm, read_num / sum(read_num_list))
    }
  }
  else if(bias == 'custom'){
    load(model)
    prob_list = c()
    for(x in transcript_polyA_arr){
      prob_list = c(prob_list, model$prob[x+1])
    }
    prob_list_norm = c()
    for(p in prob_list){
      prob_list_norm = c(prob_list_norm, p / sum(prob_list))
    }
  }
  return (prob_list_norm)
}

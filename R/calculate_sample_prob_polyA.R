#' given a transcript and fragment length, calculate probability of generating reads along the transcript
#'
#' @param transcript transcript sequence, DNAString.
#' @param fraglen fragment length, integer.
#' @return probability of generating reads along the transcript
#' @export
calculate_sample_prob_polyA = function(transcript, fraglen, polyAnum=15) {
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
  polyA_model_arr = c()

  # read in trained model:  distance to polyA, num of reads
  file_path = system.file("extdata", "summary_data_polyA_15.txt", package="Simulator")
  for(line in readLines(file_path)){
    polyA_model_arr = c(polyA_model_arr, as.numeric(line))
  }
  read_num_list = c()
  for(x in transcript_polyA_arr){
    # python index -> r index
    read_num_list = c(read_num_list, polyA_model_arr[x+1])
  }
  prob_list_norm = c()
  for(read_num in read_num_list){
    prob_list_norm = c(prob_list_norm, read_num / sum(read_num_list))
  }
  return (prob_list_norm)
}

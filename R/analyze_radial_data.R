#' Analyze output from my CellProfiler micropattern project
#'
#' @param input_file The csv output file from my micropattern project.
#' @param DNA_channel Specify the name of the DNA channel if different from DNA.
#'
#' @export
#'
#'

analyze_radial_data <- function(input_file, DNA_channel="DNA"){
  data <- utils::read.table(input_file, header=TRUE, sep=",")
  nsamples <- dim(data)[1]
  relevant_data <- as.data.frame(t(data[,stringr::str_detect(colnames(data), "MeanFrac")]))
  relevant_data$target <- stringr::str_extract(rownames(relevant_data), "(?<=MeanFrac_)[[:alnum:]]+")
  relevant_data$position <- as.numeric(stringr::str_extract(rownames(relevant_data), "(?<=_)[[:digit:]]+(?=of)"))

  max_position <- stringr::str_extract(rownames(relevant_data)[1], "(?<=of)[[:digit:]]+")
  ##data_long <- pivot_longer

  data_bg_subtracted <- relevant_data
  for (target in unique(data_bg_subtracted)){
    for (sample in 1:nsample){
      data_bg_subtracted[data_bg_subtracted$target == target, sample] <- data_bg_subtracted[data_bg_subtracted$target==target, sample] - data_bg_subtracted[data_bg_subtracted$target==target & data_bg_subtracted$position==max_position, sample]
    }
  }

  data_normalized <- data_bg_subtracted[data_bg_subtracted$position < (max_position-1),]
  for (sample in 1:nsample){
    for (position in 1:(max_position-2)){
      data_normalized[data_normalized$position==position, sample] <- data_normalized[data_normalized$position==position, sample] / data_normalized[data_normalized$position==position & data_normalized$target==DNA_channel, sample]
    }
  }
  data_normalized <- data_normalized[data_normalized$target != DNA_channel,]
  return(data_normalized)
}

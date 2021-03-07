

#' evaluate the similarities between the two set of DE genes
#'
#' @param real real data
#' @param simulated simulated data
#'
#' @return confusion matrix, logFC of the set of genes from real and simiulated data
#' @export
eval_signal  <- function(real , sim ) {


  # subset the dataset to two largest cell types
  temp_result <- process_data(real,sim )

  # find the amount of biological signals between two cell types
  result <- evaluate_DE(temp_result$real, temp_result$sim)

  return (result)
}

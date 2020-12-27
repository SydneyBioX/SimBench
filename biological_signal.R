
source("generate_DE.R")
source("evaluate_DE.R")


#' evaluate the similarities between the two set of DE genes 
#'
#' @param real real data
#' @param simulated simulated data 
#'  
#' @return confusion matrix, logFC of the set of genes from real and simiulated data
#'
eval_signal  <- function(real , sim ) {
  result_final <- list()
  
  # subset the dataset to two largest cell types 
  temp_result <- process_data(real,sim )
  
  
  # find the amount of differential expression between two cell types 
  marker <- get_DE( temp_result$real_final, temp_result$simulated_final )
  
  # genes with adjusted P.value < 0.1 are considered as DE genes 
  temp_result <- p_fdr(  real = marker$real_markers,  simulated = marker$sim_markers ,  original = real )
  
  # get the confusion matrix 
  result <- eval_metric(real_DE = temp_result$real,
                        simulated_DE =  temp_result$simulated ,
                        original =  temp_result$original  )
  
  result_final[[ "confusion" ]]  <- result
  result_final[[ "real_DE" ]]  <- temp_result$real[, c(1,2,3,4,6)]
  result_final[[ "sim_DE" ]]  <- temp_result$simulated[, c(1,2,3,4,6)]
  
  return (result_final)
}

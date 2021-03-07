

#' Evaluate DE
#'
#' @param real  A list containing real and simulated data
#' @param simulated number of cores for parallel computing
#'
#' @return DE
#' @export
evaluate_DE <- function(real, simulated) {
  # generate DE genes


  types <- c("DE", "DV", "DD", "DP", "BD")

  df <- NULL

  for (thistype in types){

    simulated_prop <- 0
    real_prop <- 0

    try({
        real_prop  <-  generate_DE( real, real$celltype, thistype)
        temp <- data.frame(  types = thistype, prop = real_prop, sim = "real")
        df <- rbind(df, temp)
    })

    try({
      simulated_prop <-  generate_DE( simulated,  simulated$celltype, thistype)
      temp <- data.frame(  types = thistype, prop =  simulated_prop, sim = "sim")
      df <- rbind(df, temp)
    })
  }

  return (df)

}








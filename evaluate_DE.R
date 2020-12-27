


library(caret)


#' evaluation similarities between the set of DE genes between real and simulated dataset
#'
#' @param real_DE the set of DE and non-DE genes from real data
#' @param simulated_DE the set of DE and non-DE genes from simulated data
#' @param original the real data used for reference  
#' 
#' @return evaluation result 
eval_metric  <- function(real_DE, simulated_DE, original ) {
  
  DE_real <-  do.call("c",  mclapply(rownames( original ) , function(x){  # loop through each gene 
    if (x %in% rownames( real_DE )  ){         
      if (real_DE[x, ]$fdr < 0.1){   
        1     # if this gene has FDR < 0.1 ,  then consider it as true DE 
      }else{
        0      # it FDR >  0.1, then consider is as true non- DE 
      }
    }else{
      0     #if Seurat did not even include this gene, it means it must be none DE 
    }
    
  }, mc.cores=4) )
  
  
  # repeat the same procedure for simulated dataset 
  DE_simulated  <-  do.call("c",  mclapply(rownames( original ) , function(x){
    if (x %in% rownames( simulated_DE)  ){
      if (simulated_DE[x, ]$fdr < 0.1){
        1
      }else{
        0
      }
    }else{
      0
    }
    
  }, mc.cores=4) )
  
  
  # provide predicted class, provide true class , positive result is the true DE 
  confused <- confusionMatrix( as.factor(DE_simulated  ) , as.factor(DE_real) , positive = "1" )
  confused <- data.frame( score = confused$byClass )
  
  
  return(confused)
}


#' perform false discovery rate p-value adjustment 
#'
#' @param real_DE the set of DE and non-DE genes from real data
#' @param simulated_DE the set of DE and non-DE genes from simulated data
#' @param original the real data used for reference  
#' 
#' @return p-value adjusted genes
p_fdr <- function ( real , simulated, original ) {
  
  simulated$fdr <- p.adjust(simulated$p_val, method= "fdr")
  real$fdr <-  p.adjust(real$p_val, method= "fdr")

  return_list <- list(real =  real, simulated = simulated, original = original)
  
  return (return_list )
  
}



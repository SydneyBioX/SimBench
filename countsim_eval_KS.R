
library(dplyr)
source( "calculateMeanVarLibrary.R" )
source( "calculateGeneCellCorr.R" )
source( "KDE_test.R" )


#' Calculate KDE test statistic on each parameter 
#'
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing 
#' 
#' @return the KDE test statistic across 13 parameters, and the raw values that are used to calculate the KDE test statistics 
countsim_eval <- function( sim_list, ncore = 8  ) {
  
  print("finding mean and variance ")
  
  obj <- calculateMeanVarLibrary(sim_list =  sim_list ,  ncore = ncore )
  
  
  
  print("calculating correlation estimate")
  maxNForCorr = 500  # default value is 500 
  # if there are more than 500 cells in this cell type, just calculate using 500 cells 
  #ie. sample - sample correlation for top 500 variable genes, across 500 samples
  # then gene - gene correlation for top 500 variable genes , across  all samples
  
  print("cellwise correlation estimate")
  sampleCorrDF <-  calculateSampleCorrs(sim_list = obj,  maxNForCorr  = maxNForCorr,  ncore = ncore)
  
  print("genewise correlation estimate  ")
  ## -------------- Calculate between-feature correlations --------------- ##
  featureCorrDF <-  calculateFeatureCorrs(sim_list = obj, ncore = ncore)
  
  
  sampleDF <- lapply(obj, function(x) {
    
    if (! is.null( x$dge$counts ) ){
      data.frame(
        Libsize = colSums(x$dge$counts),       
        Fraczero = colMeans(x$dge$counts == 0),
        TMM = x$dge$samples$norm.factors,
        type = "raw"
      ) %>% dplyr::mutate(EffLibsize = Libsize * TMM)
    }else{
      data.frame(
        Libsize = -100 , 
        Fraczero = colMeans(x$dge$log2cpm == 0) ,
        TMM =  -100 ,
        type = "normalized", 
        EffLibsize = -100        # if normalised, then can't calculate libsize or TMM
      )  
    }
  } )
  
  ns <- sapply(sampleDF, nrow)
  sampleDF <- do.call(rbind, sampleDF) %>%   dplyr::mutate(dataset = rep(names(sampleDF), ns))
  
  
  
  featureDF <- lapply(obj, function(x) {
    
    if (! is.null( x$dge$counts ) ){
      data.frame(
        average_log2_cpm = x$dge$mean,
        variance_log2_cpm = x$dge$dispersion, 
        variance_scaled_log2_cpm = x$dge$dispersion.scaled, 
        Fraczero = rowMeans(x$dge$counts == 0), 
        type = "raw"
      )
    } else {
      data.frame(
        Fraczero = rowMeans( x$dge$log2cpm == 0 ), 
        average_log2_cpm = x$dge$mean,
        variance_log2_cpm = x$dge$dispersion, 
        variance_scaled_log2_cpm = x$dge$dispersion.scaled, 
        type = "normalized"
      )
    }
  } )
  
  ns <- sapply(featureDF, nrow)
  featureDF <- do.call(rbind, featureDF) %>%   dplyr::mutate(dataset = rep(names(featureDF), ns))
  
  
  
  
  ## ----------------- Summarize data set characteristics ---------------- ##
  datasetDF <- do.call(rbind, lapply(obj, function(x) {
    
    if (! is.null( x$dge$counts ) ){ 
      data.frame(
        nVars = nrow(x$dge$counts),
        nSamples = ncol(x$dge$counts),
        type = "raw"
      ) 
      
    } else {
      data.frame(
        nVars = nrow(x$dge$log2cpm),
        nSamples = ncol(x$dge$log2cpm),
        type = "normalized"
      ) 
      
    }
  })) %>%  dplyr::mutate(dataset = names(obj))
  
  
  subsampleSize  = 10000 # choose 10000 values to calculate the similarities between the real and the simulated   
  
  dataset <- unique(  featureDF$dataset)
  real_dataset <- dataset[1]
  sim_dataset <- dataset[ 2]
  type <-   datasetDF$type[2]
  
  
  
  tmm <- efflibsize <- libsize <-  NULL
  average_log2_cpm <- variance_log2_cpm  <-  variance_scaled_log2_cpm <- NULL
  samplecor <- featurecor <- NULL
  fraczerogene <-  fraczerocell <-  NULL
  mean_variance <-   mean_fraczero <-  libsize_fraczero  <- NULL
  
  
  # can only calculate library size estimates for raw count 
  if ( type == "raw") {
   
    libsize <-  KDE_test(df =   sampleDF, column = "Libsize", subsampleSize = subsampleSize )
    
    libsize_fraczero <- KDE_test(df = sampleDF , column = c("Libsize", "Fraczero"),
                               subsampleSize = subsampleSize )
    
    tmm <-  KDE_test(df =  sampleDF, column = "TMM",    subsampleSize = subsampleSize )
    
    efflibsize <-  KDE_test(df =  sampleDF, column = "EffLibsize", 
                          subsampleSize = subsampleSize )
    
  } 
  
  # now start evaluting the stats that can be evaluated by all  
  
 
  mean_variance <-  KDE_test(df =  featureDF, column = c("average_log2_cpm", "variance_log2_cpm"), 
                           subsampleSize = subsampleSize )
  
  
  mean_fraczero <-  KDE_test(df =  featureDF , column = c("average_log2_cpm", "Fraczero"), 
                           subsampleSize = subsampleSize )
  
  
  fraczerogene <- KDE_test(df = featureDF, column = "Fraczero", 
                         subsampleSize = subsampleSize )
  
  
  average_log2_cpm <-  KDE_test(df =featureDF, column = "average_log2_cpm", 
                              subsampleSize = subsampleSize )
  
  variance_log2_cpm <-  KDE_test(df = featureDF, column = "variance_log2_cpm", 
                               subsampleSize = subsampleSize )
  
  variance_scaled_log2_cpm  <- KDE_test(df = featureDF, column = "variance_scaled_log2_cpm", 
                                      subsampleSize = subsampleSize )
  
  
 
  fraczerocell <- KDE_test(df =  sampleDF, column = "Fraczero", 
                         subsampleSize = subsampleSize )
  
  
  samplecor <- KDE_test(df =   sampleCorrDF, column = "Correlation",
                      subsampleSize = subsampleSize )
  
 
  featurecor <-  KDE_test(df =  featureCorrDF, column = "Correlation",
                        subsampleSize = subsampleSize )
  
  
  
  result  <- list()
  
  
  # can only be evaluated on raw count 
  result$tmm <-  tmm
  result$efflibsize <- efflibsize
  result$libsize <- libsize
  result$libsize_fraczero <- libsize_fraczero
  
  
  # can be evaluted using log2cpm 
  result$average_log2_cpm <-  average_log2_cpm
  result$variance_log2_cpm <-  variance_log2_cpm
  result$variance_scaled_log2_cpm <- variance_scaled_log2_cpm
  
  result$samplecor <- samplecor
  result$featurecor <-  featurecor
  
  result$fraczerogene <- fraczerogene
  result$fraczerocell <- fraczerocell 
  
  result$mean_variance <-  mean_variance
  result$mean_fraczero <- mean_fraczero
  
  
  
  result <- list(score = result, 
                 raw_value = list ( sampleDF = sampleDF,
                                    featureDF = featureDF,
                                    sampleCorrDF = sampleCorrDF, 
                                    featureCorrDF = featureCorrDF))
  return(   result)
  
}









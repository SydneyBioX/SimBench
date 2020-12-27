#' cellwise correlation 
#'
#' @param sim_list  a list containing the real and simulated data
#' @param maxNForCorr number of cells selected to compute the cellwise correlation  
#' @param ncore number of cores for parallel computing
#'
#' @return A data frame with cellwise correlations for real and simulated data
#'
calculateSampleCorrs <- function(sim_list , maxNForCorr = 1000, ncore = 8 ) {
  
  sampleCorrDF <- mclapply(sim_list , function(x) {

    # get the log2 CPM
    cpms <- as.matrix(x$dge$log2cpm)
    
    # get the top 500 most variable genes 
    var_genes <- apply(cpms, 1, var)
    maxgene <- min(500,  nrow(cpms))
    select_var <-  names(sort(var_genes, decreasing = TRUE))[1:maxgene]
    
    # Subset logcounts matrix
    cpms <- cpms[select_var, ]
    
    # say, there are 2000 cells, then only take the 1000 cells for the evaluation
    if (ncol(cpms) > maxNForCorr) {
      cpms <-  cpms[, sample(seq_len(ncol(cpms)), maxNForCorr, replace = FALSE)]
    }
    
    # although this can compute correlation using multiple cores
    # the main function itself is already set up in parallel
    # so no more multiple cores here
    corrs <- WGCNA::cor(cpms,  use =  "pairwise.complete.obs",   
                        method = "spearman",  nThreads = 1)
    
    data.frame(Correlation = corrs[upper.tri(corrs)])
    
  } ,  mc.cores = ncore )
  
  ## Merge correlations from all data sets
  ns <- vapply(sampleCorrDF, nrow, 0)
  do.call(rbind, sampleCorrDF) %>% dplyr::mutate(dataset = rep(names(sampleCorrDF), ns))
  
}




#' genewise correlation 
#'
#' @param sim_list  a list containing the real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return A data frame with genewise correlations for real and simulated data
#'
calculateFeatureCorrs <- function(sim_list ,  ncore = 8) {
  
    featureCorrDF <- mclapply(sim_list, function(x) {
      
      # get the log2 CPM
      cpms <- as.matrix(x$dge$log2cpm)
      cpms <- cpms[genefilter::rowVars(cpms) > 0,]
 
      # subset to calculate the correlation of the highly variable genes
      var_genes <- apply(cpms, 1, var)
      maxgene <- min(500,  nrow(cpms))
      select_var <-  names(sort(var_genes, decreasing = TRUE))[1:maxgene]
      cpms <- cpms[select_var, ]
    
      corrs <- WGCNA::cor(t(cpms),   use =  "pairwise.complete.obs",
                          method = "spearman" ,  nThreads =  1)
      
      data.frame(Correlation = corrs[upper.tri(corrs)])

    } ,  mc.cores = ncore)
    
    ## Merge correlations from all data sets
    ns <- vapply(featureCorrDF, nrow, 0)
    do.call(rbind, featureCorrDF) %>%  dplyr::mutate(dataset = rep(names(featureCorrDF), ns))
    
  }




#' Process data
#'
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return the KDE test statistic across 13 parameters, and the raw values that are used to calculate the KDE test statistics
#' @export
process_data <- function(real, simulated)   {


  total_celltypes <- unique ( as.character( real$celltype ))
  if (length(total_celltypes) < 2){
    stop("not enough cell types")
  }

  # get the two most abundant cell types to evaluate
  celltype <- names( sort(table(  real$celltype),decreasing=TRUE)[1:2]  )

  # subset to the two most abundant cell types
  real_thiscelltype <- real[ , real$celltype  %in%  celltype  ]
  sim_thiscelltype <-  simulated[ ,  simulated$celltype %in% celltype ]


  # if the dataset is singlecellexperiment, change it to seurat
  # this is because we will get the DE genes using Seurat
  if (class(real_thiscelltype) == "SingleCellExperiment" ) {
    temp_celltype <- real_thiscelltype$celltype
    real_thiscelltype <- CreateSeuratObject( counts = counts(real_thiscelltype))
    real_thiscelltype$celltype <- temp_celltype
    real_thiscelltype <- NormalizeData(real_thiscelltype)
  }

  if (class( sim_thiscelltype) == "SingleCellExperiment"  ){
    temp_celltype <- sim_thiscelltype$celltype
    sim_thiscelltype <- CreateSeuratObject( counts = counts(sim_thiscelltype))
    sim_thiscelltype$celltype <- temp_celltype
    sim_thiscelltype <- NormalizeData(sim_thiscelltype)
  }



  return (list(real_final = real_thiscelltype  ,  simulated_final = sim_thiscelltype)  )
}



#' Generate DE
#'
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return the KDE test statistic across 13 parameters, and the raw values that are used to calculate the KDE test statistics
#' @export
generate_DE <- function(exprsMat,
                             trainClass,
                             feature = c("DE", "DV", "DD", "DP", "BD"),
                             # limma is DE ,  chisq is DP , BI is Bimodal distribution
                             pSig = 0.1
){


  if (feature == "DV") {
    tt <- doDV(exprsMat, trainClass)
    tt <-   sum(tt < pSig)

  } else if (feature == "DD") {
    tt <- doDD(exprsMat, trainClass)
    tt <-   sum(tt < pSig)

  } else if (feature == "DP") {
    tt <- doChisSquared(exprsMat, trainClass)
    tt <-   sum(tt < pSig , na.rm = T)

  } else if (feature == "BD") {
    tt <- doBI(exprsMat, trainClass)
    tt <-   sum(tt > 0.3 , na.rm = T)

  }

  else{
    tt <- doLimma(exprsMat, trainClass)
    tt <-  sum(tt$adj.P.Val < pSig)
  }
  res <- tt / dim( exprsMat )[1]

  return(res)
}


#' @importFrom limma eBayes lmFit
#' @importFrom methods new
doLimma <- function(exprsMat, cellTypes, exprs_pct = 0.05){

  cellTypes <- droplevels(as.factor(cellTypes))


  tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[1], 1, 0))
  design <- stats::model.matrix(~tmp_celltype)


  meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
    Matrix::rowMeans(exprsMat@assays$RNA@data[, tmp_celltype == i, drop = FALSE])
  }))

  meanPct <- do.call(cbind, lapply(c(0,1), function(i){
    Matrix::rowSums(exprsMat@assays$RNA@data[, tmp_celltype == i,
                                             drop = FALSE] > 0)/sum(tmp_celltype == i)
  }))

  keep <- meanPct[,2] > exprs_pct

  y <- methods::new("EList")
  y$E <- exprsMat@assays$RNA@data[keep, ]
  fit <- limma::lmFit(y, design = design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt   <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)



  if (!is.null(tt$ID)) {
    tt <- tt[!duplicated(tt$ID),]
    rownames(tt) <- tt$ID
  }

  tt$meanExprs.1 <- meanExprs[rownames(tt), 1]
  tt$meanExprs.2 <- meanExprs[rownames(tt), 2]
  tt$meanPct.1 <- meanPct[rownames(tt), 1]
  tt$meanPct.2 <- meanPct[rownames(tt), 2]


  return(tt)


}


doDV <- function(exprsMat, cellTypes){


  cellTypes <- droplevels(as.factor(cellTypes))

  tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[1], 1, 0))

  meanPct <- do.call(cbind, lapply(c(0,1), function(i){
    Matrix::rowSums(exprsMat@assays$RNA@data[,
                                             tmp_celltype == i,
                                             drop = FALSE] > 0)/sum(tmp_celltype == i)
  }))


  posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
  # print(sum(posNeg))
  exprsMat_filt <- exprsMat@assays$RNA@data[posNeg,]
  tt <- apply(exprsMat_filt, 1, function(x) {
    df <- data.frame(gene = x, cellTypes = as.factor(tmp_celltype))
    stats::bartlett.test(gene~cellTypes, df)$p.value
  })

  tt  <- stats::p.adjust(tt , method = "BH")




  return(tt)


}




doDD <- function(exprsMat, cellTypes){

  cellTypes <- droplevels(as.factor(cellTypes))


  tmp_celltype <- ifelse(cellTypes == levels(cellTypes)[1], 1, 0)

  meanPct <- do.call(cbind, lapply(c(0,1), function(i){
    Matrix::rowSums(exprsMat@assays$RNA@data [,
                                              tmp_celltype == i,
                                              drop = FALSE] > 0)/sum(tmp_celltype == i)
  }))

  posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05

  exprsMat_filt <- exprsMat@assays$RNA@data[posNeg,]

  tt  <- apply(exprsMat_filt, 1, function(x) {
    x1 <- x[tmp_celltype == 0]
    x2 <- x[tmp_celltype == 1]
    stats::ks.test(x1, x2, alternative = "greater")$p.value
  })



  tt <- stats::p.adjust(tt , method = "BH")

  return(tt)


}




doChisSquared <- function(exprsMat, cellTypes, threshold = 1){


  cellTypes <- droplevels(as.factor(cellTypes))


  tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[1], 1, 0))


  zerosMat <- ifelse(as.matrix( exprsMat@assays$RNA@data  ) > threshold, 1, 0)

  tt  <- apply(zerosMat,1,  function(x){
    tab <- c()
    for (i in c(0,1)) {
      tmp <- factor(x[tmp_celltype == i], levels = c(0, 1))
      tab <- rbind(tab, table(tmp))
    }

    suppressWarnings(stats::chisq.test(tab)$p.value)

  })

  tt <- stats::p.adjust(tt , method = "BH")

  return(tt)


}


doBI <- function(exprsMat, cellTypes){
  # Select genes by bimodal index

  cellTypes <- droplevels(as.factor(cellTypes))


  tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[1], 1, 0))

  pi <- table(tmp_celltype)/length(tmp_celltype)

  agg_mean <- do.call(cbind, lapply(c(0,1), function(i){
    Matrix::rowMeans(exprsMat@assays$RNA@data[, tmp_celltype == i, drop = FALSE])
  }))

  agg_sd2 <- do.call(cbind, lapply(c(0,1), function(i){
    apply(exprsMat@assays$RNA@data[, tmp_celltype == i, drop = FALSE], 1, stats::var)
  }))

  bi <- abs(agg_mean[,2] - agg_mean[,1])/sqrt(pi[1]*agg_sd2[,1] +
                                                pi[2]*agg_sd2[,2])

  bi <- unlist(bi)
  names(bi) <- rownames(exprsMat)
  bi <- bi[order(bi, decreasing = TRUE)]
  tt  <- bi


  return(tt)


}

#' Calculate mean variance and library size estimates
#'
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return a list containing the mean , variance , library size estimates of the real and simulated data for downstream evaluation
#' @import Seurat
#' @import edgeR
#' @import DESeq2
#' @export
calculateMeanVarLibrary <- function(sim_list , ncore = 8) {

  mclapply(sim_list, function(ds) {
    # if Seurat , it means the simulated count matrix is normalised
    # hence library size estimate cannot be obtained
    if (class(ds) == "Seurat") {
      dge <- list()
      dds <- list()

      msg <- try({
        temp_seurat <-  FindVariableFeatures(ds, selection.method = "disp")
      })

      if (class(msg) == "try-error") {
        ds <- NormalizeData(ds, scale.factor = 1e6)
      }

      temp_seurat <- FindVariableFeatures(ds, selection.method = "disp")
      dge$log2cpm <- ds[["RNA"]]$data
      dge$mean  <-  HVFInfo(temp_seurat)[, 1]
      dge$dispersion <-   HVFInfo(temp_seurat)[, 2]
      dge$dispersion.scaled <- HVFInfo(temp_seurat)[, 3]


    } else if (class(ds) == "DESeqDataSet") {
      # if the object is DESeq, it means we can find the library size estimate

      dge <- edgeR::DGEList(counts = DESeq2::counts(ds))
      ## Calculate normalization factors
      dge <- edgeR::calcNormFactors(dge)


      # if the class is DEseq, it means the data is raw, and need cpm normalization
      temp_seurat <-  CreateSeuratObject(counts = counts(ds))
      temp_seurat <-  NormalizeData(temp_seurat,  scale = 1e6)
      temp_seurat <- SetAssayData(object =   temp_seurat, slot = "counts", new.data =  temp_seurat[["RNA"]]$data)

      # note that even though it "finds variable features"
      # it returns the mean and sd of all genes in the dataset
      temp_seurat <-  FindVariableFeatures(temp_seurat, selection.method = "disp")
      dge$mean  <-  HVFInfo(temp_seurat)[, 1]
      dge$dispersion <-  HVFInfo(temp_seurat)[, 2]
      dge$dispersion.scaled <-  HVFInfo(temp_seurat)[, 3]
      dge$log2cpm <- temp_seurat[["RNA"]]$data



      ## --------------------------- DESeq2 -------------------------- ##
      ## Calculate size factors
      dds <- DESeq2::estimateSizeFactors(ds, type = "poscounts")
    }

    # the dge is used to store mean and variance
    # dds is used to store library size
    list(dge = dge, dds = dds)

  } , mc.cores = ncore)

}

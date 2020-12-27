

library(Seurat)
library(SingleCellExperiment)

#' Process the data into Seurat object and select the two most abundant cell types
#'
#' @param real real data
#' @param simulated simulated data 
#'  
#' @return processed data
#'
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






#' obtain log fold change of each gene in the real and simulated data
#'
#' @param real real data
#' @param simulated simulated data 
#'  
#' @return log fold change of each gene
#'
get_DE <- function(real, simulated){
  
  celltype <- unique( real$celltype )
  
  
  real_markers <- FindMarkers(real, ident.1 =   celltype[1], group.by = 'celltype' ,
                              logfc.threshold = 0,  test.use = "wilcox", min.pct = 0)
  
  sim_markers <- FindMarkers(simulated, ident.1 =   celltype[1], group.by = 'celltype' ,
                             logfc.threshold = 0,  test.use = "wilcox", min.pct = 0)
  
  
  return (list( real_markers =  real_markers ,  sim_markers  =  sim_markers )  ) 
}













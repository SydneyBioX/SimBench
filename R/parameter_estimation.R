

#' format the data for evaluation
#'
#' @param sim_list  a list for storing the real and simulated data
#' @param thissim the data to be put into the sim_list
#' @param name name of the data (eg, "real", "simulated")
#' @param celltype cell type that the data contains
#' @param counttype whether the data is unnormalised or normalised
#' @import DESeq2
#' @return a sim_list filled with the formatted data
#' @export
prepare_eval <- function(  sim_list, thissim , name , celltype,  counttype ){


  # if count type is normalized , then put it in Seurat
  # because DESeqDataSetFromMatrix cannot handle normalized count with decimals
  if ( counttype == "normalized"  ){

    if (class(thissim) == "SingleCellExperiment"){
      meta_data <-  as.data.frame( colData(thissim))
      meta_data <- meta_data[, -1]
      thissim <-  CreateSeuratObject(counts = counts(thissim ),   meta.data =    meta_data)
      thissim$group <- celltype
      thissim$sample <- colnames( thissim )
    }

  } else if (counttype == "raw"){

    # if count type is raw , then put it in a DESeqDataSetFromMatrix
    count <-  as.matrix(  counts( thissim)  )
    mode( count) <- "integer"

    thissim  <- DESeqDataSetFromMatrix(countData =  count,
                                       colData =  data.frame(  group = celltype ,
                                                               sample =  colnames(thissim )   ,
                                                               row.names = colnames(thissim )  ,
                                                               stringsAsFactors = FALSE) ,
                                       design = ~ 1)

  }
  sim_list[[name ]] <- thissim


  return(sim_list)
}






#' summarise the scores across multiple cell type according to the proportion of each cell type
#'
#' @param intermediate  the scores from multiple cell type
#'
#' @return a summarised score
#' @export
summarise_score <- function( intermediate   ){

  return_dataframe <- NULL

  allparameter <- names(intermediate)

  for ( i in (1:length(allparameter)) ){

    thisparameter <- intermediate[[i]]
    thisparameter$temp <- as.numeric( thisparameter$kde_zstat ) * thisparameter$proportion

    temp <-  thisparameter %>%  dplyr::summarise( sum_kde = sum( temp)) 

    return_dataframe[[  allparameter[i]]] <- temp
  }

  return (  return_dataframe )

}






#' run evaluation on parameter estimation
#'
#' @param real real data
#' @param sim simulated data
#' @param type  whether the simulated data is raw (unnormalised raw count) or normalized
#' @param method name of the simulation method

#' @return distribution of parameters that are used to construct KDE test, KDE test statistics of each cell type, combined KDE test statistics
#' @export

eval_parameter <- function(real  , sim  , type  = "raw"  ,  method = "samplemethod") {

    #-----------Prepare the dataset ------------------#

    celltype <- unique( as.character( real$celltype))

    eval_result_allcelltype <- NULL
    eval_result_raw <- list()


    #----------------- start evaluating --------------#

    # loop through each cell type
    for (i  in  (1: length( celltype) )) {
      sim_list <- list()

      thiscelltype <-  celltype[i]
      celltypeproportion <- sum( real$celltype == thiscelltype ) / length( real$celltype)
      print(paste0 ( "evaluating cell type .. " ,  thiscelltype   ))

      # prepare the format
      data_real <-  real[ ,  real$celltype ==   thiscelltype ]
      if(ncol(data_real) <10 ) next
      sim_list <-  prepare_eval( sim_list = sim_list,
                                 thissim =  data_real ,
                                 name =   "real" ,
                                 celltype = thiscelltype,
                                 counttype =  type)

      data_sim  <-  sim [ , sim$celltype ==   thiscelltype ]

      sim_list <-  prepare_eval( sim_list  = sim_list,
                                 thissim  = data_sim ,
                                 name = method  ,
                                 celltype = thiscelltype,
                                 counttype = type )


      # generate the evaluation result
      eval_result <- countsim_eval( sim_list )

      eval_result_raw[[thiscelltype]]  <-  eval_result

      eval_result <- eval_result$score


      # add additional column on cell type information
      eval_result  <- lapply( eval_result  , function(x)
        cbind(x,  celltype = thiscelltype))

      # add additional column on  cell type proportion information
      eval_result  <- lapply( eval_result  , function(x)
        cbind(x,  proportion =     celltypeproportion))


      # combine the result for each cell type
      if   ( is.null( eval_result_allcelltype )  ){
        eval_result_allcelltype <-  eval_result
      } else{
        eval_result_allcelltype <-  mapply(FUN = rbind, eval_result_allcelltype ,   eval_result, SIMPLIFY = FALSE)
      }

    } #end the for loop for evaluating each individual cell type


    summarise_celltype <- summarise_score(eval_result_allcelltype)

    return (  list ( raw_value =  eval_result_raw ,  stats_celltype = eval_result_allcelltype ,
                    stats_overall =  summarise_celltype ) )

}

# save the result for the methods that evaluate based on each cell type


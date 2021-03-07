

#' Perform KDE test
#'
#' @param df dataframe containing the raw distribution for comparison
#' @param column  the column containing the parameter of interest
#' @param subsampleSize maximum number of data points for comparison
#'
#' @return KDE test statistic
#' @import ks
#'
#' @export
KDE_test <- function(df,  column,  subsampleSize ) {


    ## Initialize data frame by populating it with all data set pairs
      method_name <-  unique(df$dataset)

      ds1 =  method_name[1]    #set the name, eg, "real"
      ds2 =  method_name[2]   # set the name, eg, "splatter"


      ## Remove rows with NA values in column(s) of interest
      df <- df[rowSums(is.na(df[, column, drop = FALSE])) == 0, ]

      # kde
      if (nrow(df)  >  subsampleSize) {
        df <-   df[sample(seq_len(nrow(df)),  subsampleSize, replace = FALSE), ]
      }

      # KDE test requires vector for univariate and matrix for multivariate comparison
      if (length(column) == 1) {
        kde <- kde.test(x1 = as.vector(df[, column][df$dataset == ds1]) ,
                        x2 = as.vector(df[, column][df$dataset == ds2]))
      } else{
        kde <- kde.test(x1 = as.matrix(df[, column][df$dataset == ds1 , ]) ,
                        x2 = as.matrix(df[, column][df$dataset == ds2, ]))
      }

      ds_res <- data.frame( kde_tstat =  kde$Tstat ,
                               kde_zstat = kde$zstat ,
                               kde_pvalue = kde$pvalue)

    ## Return output table
    return( ds_res)

}

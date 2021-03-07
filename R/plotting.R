#  the key challenge is to compare across each evaluation criteria.
# For this, show the distribution of different metrics


# show that for the same dataset, the score share similar distribution, ie, needs to be comparable



#' Plot
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return the KDE test statistic across 13 parameters, and the raw values that are used to calculate the KDE test statistics
#' @import ggplot2 ggpubr ggthemes
#' @export
draw_parameter_plot <- function( result  ){

      sampleDF <- result$sampleDF
      featureDF <- result$featureDF
      sampleCorrDF <- result$sampleCorrDF
      featureCorrDF <- result$featureCorrDF


      plot_list <- list()

      th <-   theme(text=element_text(size=12 ),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=0.2, fill=NA) )


        p <-  ggplot( sampleDF , aes(x = Libsize , group = dataset,  fill=dataset , color = dataset )) +
          geom_density( alpha = 0.7 ) +
          xlab("library size")  +
          scale_fill_manual(values=c( "#184275", "#b3202c" )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle( "libsize") + th

        plot_list$libsize <- p


        p <-  ggplot( sampleDF , aes(x =  Libsize , y = Fraczero , color = dataset )) +
          geom_point(size = 0.5, alpha = 0.5 ) +
          xlab("library size") + ylab("fraction zero per gene")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c" )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("libsize_fraczero")+ th

        plot_list$libsize_fraczero <- p


        p <-  ggplot( sampleDF , aes(x = TMM , group = dataset,  fill=dataset, color = dataset  )) +
          geom_density( alpha = 0.7 ) +
          xlab("TMM")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("TMM") + th

        plot_list$tmm <- p



        p <-  ggplot( sampleDF , aes(x = EffLibsize, group = dataset,   fill=dataset, color = dataset  )) +
          geom_density( alpha = 0.7 ) +
          xlab("effective library size")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("effective library size") + th

        plot_list$effectivelibsize <- p




        p <-  ggplot( featureDF  , aes(x =  average_log2_cpm , y = variance_log2_cpm , color = dataset, fill=dataset )) +
          geom_point(size = 0.5, alpha = 0.1) +
          xlab(" mean expression ") + ylab(" variance of gene expression ")    +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle( "mean_variance" ) + th

        plot_list$mean_variance  <- p


        p <-  ggplot(featureDF, aes(x = variance_log2_cpm , group = dataset,   fill=dataset , color = dataset )) +
          geom_density( alpha = 0.7 ) +
          xlab("variance log2 cpm")    +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("variance") + th

        plot_list$variance  <- p




        p <-  ggplot(featureDF, aes(x = variance_scaled_log2_cpm  , group = dataset, fill=dataset , color = dataset )) +
          geom_density( alpha = 0.7 ) +
          xlab("variance scaled log2 cpm")   +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("scaled variance") + th

        plot_list$variance_scaled  <- p





        p <-  ggplot( sampleCorrDF , aes(x =  Correlation, group = dataset,  fill=dataset , color = dataset )) +
          geom_density( alpha = 0.7 ) +
          xlab("sample correlation")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("samplecor")  + th

          plot_list$samplecor <- p



        p <-  ggplot(featureCorrDF , aes(x =  Correlation, group = dataset,   fill=dataset, color = dataset  )) +
          geom_density( alpha = 0.7 )  +
          xlab("feature correlation")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("featurecor") + th

        plot_list$featurecor <- p


        p <-  ggplot( featureDF  , aes(x =  average_log2_cpm , y = Fraczero , color = dataset)) +
          geom_point(size = 0.5, alpha = 0.1) +
          xlab("mean expression") + ylab("fraction zero per gene")   +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("mean_fraczero") + th

      plot_list$mean_fraczero <- p




        p <-  ggplot(featureDF, aes(x = Fraczero, group = dataset,   fill=dataset, color = dataset  )) +
          geom_density( alpha = 0.7 )  +
          xlab("Fraction zeros per gene")   +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("fraczerogene") + th

      plot_list$fraczerogene <- p



        p <-  ggplot(sampleDF, aes(x = Fraczero, group = dataset,   fill=dataset  , color = dataset )) +
          geom_density( alpha = 0.7 ) +
          xlab("Fraction zeros per cell")   +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_colour_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("fraczerocell") + th

       plot_list$fraczerocell <- p




        p <-  ggplot(featureDF, aes(x = average_log2_cpm  , group = dataset,   fill=dataset, color = dataset  )) +
          geom_density( alpha = 0.7 ) +
          xlab("average log2 cpm")  +
          scale_fill_manual(values=c(  "#184275", "#b3202c"  )) +
          scale_color_manual(values=c(  "#184275", "#b3202c"  )) +
          ggtitle("mean") + th

       plot_list$mean <- p


       return( plot_list )


}



#' Plot biological plot
#'
#' @param sim_list  A list containing real and simulated data
#' @param ncore number of cores for parallel computing
#'
#' @return the KDE test statistic across 13 parameters, and the raw values that are used to calculate the KDE test statistics
#' @import ggplot2 ggpubr ggthemes
#' @export
draw_biosignal_plot <- function(result){


   result$types <- factor(result$types , c("DE" , "DV", "DD", "DP", "BD"))

   p <-  ggplot(data= result , aes(x= types, y=prop, fill=sim)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_manual(values=c("#b3202c" ,  "#184275"  )) + theme_minimal()  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=0.2, fill=NA)) + ylim(0,1)


  return(p)


}



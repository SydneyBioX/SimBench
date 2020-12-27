

# Getting started 

This github repo contains the codes for   
i) quantifying the distributional similarities between a simulated scRNA-seq and a real scRNA-seq data using KDE test (Kernel Density Based Global Two-Sample Comparison Test)   
ii) measuring the similarities between the set of differentially expressed genes in a simulated scRNA-seq and a real scRNA-seq data. 




## Dependencies 

Running the source codes requires the following R packages: 

```
ggthemes, ggpubr, ggplot2, dplyr, plyr, Seurat, SingleCellExperiment, edgeR, DESeq2, caret, ks
```




## Example usage


We have provided a 'simulated' data (`sim.rds`) and a 'real' (`real.rds`) scRNA-seq in the github folder to illustrate the usage of our codes. 




### Load scripts 

The following are the scripts containing the main functions. 

```
source("parameter_estimation.R")  
source("plotting.R")  
source("biological_signal.R")
```

### Load example data 

The files are provided in the `data` folder in this github repo .

```
real <-  readRDS("real.rds")
sim <-  readRDS("sim.rds")

```

Note both the sim and real dataset needs to be SingleCellExperiment object.  
If `celltype` is provided in the object, then the comparison will be made based on each cell type and then combined using a weighted sum (where the weight is the proportion of the cell type).  
if no `celltype` is provided, then the comparison will be made based on the entire dataset. 



### Parameter estimation 

The parameter estimation score can be obtained by : 

```
parameter_result <- eval_parameter(real = real, sim = sim, type = "count" , method = "samplemethod")
```
The output contains 3 fields:   
`stats_overall` gives the overall KDE test statistics     
`stats_celltype` gives the KDE test statistics for each cell type   
`stats_raw` gives the raw values used to perform the KDE test (eg, the mean expression of each gene) 


#### Visualise 

We can use the raw value to visualise the simulated dataset and real dataset over 13 parameters. 

```
distribution_celltype <- parameter_result$raw_value$`B cell`$raw_value #this obtain the distribution of B cell type 
fig <- draw_plot(distribution_celltype) 
ggarrange( plotlist =  fig ,  common.legend = T)

```


### Maintaining biological signatures

  
Evaluation of biological signals can be obtained by 

```
signal_result <- eval_signal( real = real, sim = sim  )
```

The output contains 3 fields:
`confusion` gives the evaluation result on precision, recall, etc   
`real_DE` gives the log fold change and adjusted P-value of each gene in the real data
`sim_DE` gives the log fold change and adjusted P-value of each gene in the simulated data


# Reference

Part of the codes are taken from the R package `countsimQC`.    

> Code: https://github.com/csoneson/countsimQC,   
Published article: Soneson, C., & Robinson, M. D. (2018). Towards unified quality verification of synthetic count data with countsimQC. Bioinformatics, 34(4), 691-692.).   

However installation of `countsimQC` is not required for running our codes. 
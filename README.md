# histoCAT_SCNA
histoCAT from Schapiro et al. 2017 used from IMCR tools for SCNA. https://www.nature.com/articles/nmeth.4391

# Introduction

HistoCAT finds statistically enriched interactions between cell phenotypes in cellular neighborhoods. The original method uses a distance to define a neighborhood, while IMCR tools enables the user to use a Delaunay graph for neighborhood definition.Pairwise interactions between celltypes are compared to a random distribution using two individual one-tailed permutation tests. The the permutation test for comparison is performed by Monte Carlo sampling. The p-values indicate either avoidance or interaction of the cell types. HistoCAT provides direction, so that celltype A to celltype B ist not the same as celltype B to celltype A.

# Installation

HistoCAT is run via the R package imcRtools (1.4.2) which can be installed from Bioconductor

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("imcRtools")
```

## Dependencies

The script requires the R packages cytomapper (1.10.1), imcRtools (1.4.2), SingleCellExperiment (1.20.0), readr (2.1.3), SpatialExperiment (1.8.0) and tidyverse (1.3.2) to be installed. 

## Source

The script was written briefly following this tutorial by the authors https://bodenmillergroup.github.io/imcRtools/articles/imcRtools.html

# Usage

Clone the github repo and create an empty "output" folder in the same folder as the github repo was cloned into. The current script takes data in form of .csv files from a "data" folder two levels up the folder you cloned the repo in. The .csv files have three columns, one with x, one with y, and one with ct annotations. 

Make sure the structure is as described or change the script to your needs.

You can run the script from the command line when navigated into the `src` folder

`Rscript histoCAT.R`

An output .csv file with all statsitics histoCAT provides will be generated in the output folder for downstream analysis.






 
 



# IMCRtools (classic and histoCAT)
histoCAT from Schapiro et al. 2017 and IMCR classic used with IMCR tools for NEP comparison. https://www.nature.com/articles/nmeth.4391

# Introduction

HistoCAT finds statistically enriched NEPs between cell phenotypes in cellular neighborhoods. The original method uses a distance to define a neighborhood, while IMCR tools enables the user to use a Delaunay graph for neighborhood definition.Pairwise NEPs between celltypes are compared to a random distribution using two individual one-tailed permutation tests. The permutation test for comparison is performed by Monte Carlo sampling. The p-values indicate either avoidance or interaction of the cell types. HistoCAT provides direction, so that celltype A to celltype B ist not the same as celltype B to celltype A.We extracted the p_lt values for downstream comparison. 

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

The script was written following this tutorial by the authors https://bodenmillergroup.github.io/imcRtools/articles/imcRtools.html

## Usage

Load all required packages. The input path, output path, neighborhood (delaunay, knn) and method (classic, histoCAT) definitions need to be set at the top of the script.
The scripts in the folder src were used to create the results on simulated data (src/histoCAT_simulation.R) and the myocardial infarction dataset from Wuennemann et al. (src/histoCAT_MI.R) in the Schiller et al. (2025) manuscript on NEP comparison. 
The script was run on an M2 MacBookPro. 
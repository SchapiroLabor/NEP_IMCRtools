# IMCRtools (classic and histoCAT)
histoCAT from Schapiro et al. 2017 and IMCR classic used with IMCR tools for NEP comparison. https://www.nature.com/articles/nmeth.4391

IMCRtools is an analysis suite for IMC data analysis. It includes a function to determine statistically enriched cell-cell interactions with a classic or histoCAT version. Both find statistically enriched NEPs between cell phenotypes in cellular neighborhoods. The original histoCAT method uses a distance to define a neighborhood, while IMCR tools enables the user to use a Delaunay graph for neighborhood definition. Pairwise NEPs between celltypes are compared to a random distribution using two individual one-tailed permutation tests. We extracted the p_lt values for downstream comparison. 

# Usage

The R package imcRtools (1.4.2) can be installed from Bioconductor

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("imcRtools")
```
The script was run on an M2 MacBookPro. 
Other dependencies

- BiocParallel (>= 1.36.0)
- lubridate (>= 1.9.3)
- forcats (>= 1.0.0)
- stringr (>= 1.5.0)
- dplyr (>= 1.1.3)
- purrr (>= 1.0.2)
- tidyr (>= 1.3.0)
- tibble (>= 3.2.1)
- ggplot2 (>= 3.4.4)
- tidyverse (>= 2.0.0)
- readr (>= 2.1.4)
- imcRtools (>= 1.8.0)
- SpatialExperiment (>= 1.12.0)
- cytomapper (>= 1.14.0)
- SingleCellExperiment (>= 1.24.0)
- SummarizedExperiment (>= 1.32.0)
- Biobase (>= 2.62.0)
- GenomicRanges (>= 1.54.1)
- GenomeInfoDb (>= 1.38.0)
- IRanges (>= 2.36.0)
- S4Vectors (>= 0.40.1)
- BiocGenerics (>= 0.48.1)
- MatrixGenerics (>= 1.14.0)
- matrixStats (>= 1.1.0)
- EBImage (>= 4.44.0)

## Source

The script was written following this tutorial by the authors https://bodenmillergroup.github.io/imcRtools/articles/imcRtools.html

## Data

### In silico tissue (IST) data
Simulated .csv data with x, y, and ct annotation columns were used. The asymmetric and symmetric in silico tissue (IST) datasets were generated as described here: https://github.com/SchapiroLabor/NEP_IST_generation. 

### Myocardial infarction (MI) data

UPDATE PATH


## Scripts

`/src`:
- `/IMCRtools_simulation.R`: This script runs IMCRtools on the simulated data with a delaunay neighborhood definition. By setting the method parameter at the top of the script to "histocat" or "classic", the method histoCAT with conditional normailzation ot the classic method with total cell type count normalization are run, respectively. 
- `/IMCRtools_MI.R`: This script runs IMCRtools on the MI data using k-nearest neighbors as neighborhood definition with k=5. By setting the method parameter at the top of the script to "histocat" or "classic", the method histoCAT with conditional normailzation ot the classic method with total cell type count normalization are run, respectively. 

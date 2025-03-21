# IMCRtools (classic and histoCAT)
histoCAT from Schapiro et al. 2017 and IMCR classic used with IMCR tools for NEP comparison. https://www.nature.com/articles/nmeth.4391

IMCRtools is an analysis suite for IMC data analysis. It includes a function to determine statistically enriched cell-cell interactions with a classic or histoCAT version. Both find statistically enriched NEPs between cell phenotypes in cellular neighborhoods. The original histoCAT method uses a distance to define a neighborhood, while IMCR tools enables the user to use a Delaunay graph for neighborhood definition. Pairwise NEPs between celltypes are compared to a random distribution using two individual one-tailed permutation tests. We extracted the sigval values for downstream comparison. 

# Usage

The R package imcRtools (1.4.2) can be installed from Bioconductor

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("imcRtools")
```
The script was run on an M2 MacBookPro. 
Other dependencies are specified in `session_info.txt`


## Source

The script was written following this tutorial by the authors https://bodenmillergroup.github.io/imcRtools/articles/imcRtools.html

## Data

### In silico tissue (IST) data
Simulated .csv data with x, y, and ct annotation columns were used. The asymmetric and symmetric in silico tissue (IST) datasets were generated as described here: https://github.com/SchapiroLabor/NEP_IST_generation. 

### Myocardial infarction (MI) data

Sequential Immunofluorescence data was accessed via Synapse (project SynID : syn51449054): https://www.synapse.org/Synapse:syn51449054. The dataframe with phenotypes was accessed within the project at:  https://www.synapse.org/Synapse:syn65487454.

## Scripts

`/src`:
- `/IMCRtools_simulation.R`: This script runs IMCRtools on the simulated data with a delaunay neighborhood definition. By setting the method parameter at the top of the script to "histocat" or "classic", the method histoCAT with conditional normailzation ot the classic method with total cell type count normalization are run, respectively. 
- `/IMCRtools_MI.R`: This script runs IMCRtools on the MI data using k-nearest neighbors as neighborhood definition with k=5. By setting the method parameter at the top of the script to "histocat" or "classic", the method histoCAT with conditional normailzation ot the classic method with total cell type count normalization are run, respectively. 

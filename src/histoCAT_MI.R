#################################################
### Script for running IMCR tools for MI data ###
### Schiller et at. 2025.                     ###
### Author: Chiara Schiller                   ###
#################################################

### This script runs IMCRtools histoCAT or classic for extracting neighbor preference by sample matrices for NEP method comparison

### load packages
library(cytomapper)
library(imcRtools)
library(SingleCellExperiment)
library(readr)
library(SpatialExperiment)
library(tidyverse)
library(BiocParallel)

### Settings
# choose from histocat or classic
method = "histocat"
# choose nbh def from knn or delaunay
nbh_type = "delaunay"
# if knn, choose k
k_number = 5
# adapt colPairName name depending on neighborhood definition
if(nbh_type == "knn")
  colPairName = "knn_interaction_graph"
if(nbh_type == "delaunay")
  colPairName = "delaunay_interaction_graph"
# data paths
data_path = "./../../../../../MI_heart_paper/data/cell_table_final.csv"
output_path = "./../../../Comparison/20250218_results_MI/histoCAT_p_lt_knn5.csv"

### Data preparation
# load and prepare data for SpatialExperiment
data = read_csv(data_path)
data = data[1:1000,]
data = data %>% filter(final_cell_type != "exclude")

# create dataframe with random values to pretend to have expression data to generate a Spatial Experiment object
df = cbind(rep(1,base::nrow(data)),rep(1,base::nrow(data)))
colnames(df) = c("M1", "M2")
metadata = as.data.frame(cbind(as.character(data$final_cell_type), as.character(data$fov), rownames(data)))
colnames(metadata) = c("ct", "img_id","cell_ID")

# create Spatial Experiment
spe = SpatialExperiment(assay=as.data.frame(t(df)), 
                        colData=metadata, 
                        spatialCoords = as.matrix(data[,(names(data) %in% c("X_centroid", "Y_centroid"))])
)

### NEP analysis
# delaunay, knn or distance neighborhood definition
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = nbh_type,
                         k = k_number,
                         coords = c("X_centroid", "Y_centroid"))

assayNames(spe) = "expr"

spe <- aggregateNeighbors(spe,
                          colPairName = colPairName,
                          aggregate_by = "metadata",
                          count_by = "ct") 

out <- countInteractions(spe,
                           group_by = "img_id",
                           label = "ct",
                           method = method,
                           colPairName = colPairName)
# Parallelize the testInteractions function and set seed
BPPARAM <- MulticoreParam(workers = parallel::detectCores() - 2, RNGseed = 123)  # Adjust the number of cores
out <- testInteractions(
                          spe,
                          group_by = "img_id",
                          label = "ct",
                          method = method,
                          colPairName = colPairName,
                          BPPARAM = BPPARAM
                        )

### create standard matrix for output comparison
# replace values that are not significant with False
out$p_lt[out$sig == FALSE] <- 0.5

data = as.data.frame(out) %>% select(from_label, to_label, p_lt, group_by)
data <- data.frame(lapply(data, function(x) ifelse(is.na(x) | is.nan(x), 0, x)))
data = as.data.frame(data)
data$key = paste(data$from_label, data$to_label, sep = "_")
data = data %>% select(-c(from_label, to_label))
data = spread(data, key = key, value = p_lt)

rownames(data) = sub(".csv", "", data$group_by)
data = data[,-1]

### save output
write.csv(data,file=output_path,row.names = TRUE)

### heatmap
library(pheatmap)
pheatmap(data)


sessionInfo()


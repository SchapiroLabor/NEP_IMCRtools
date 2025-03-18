########################################################
### Script for running IMCR tools for simulated data ###
### Schiller et at. 2025.                            ###
### Author: Chiara Schiller                          ###
########################################################

### This script runs IMCRtools histoCAT or classic for extracting neighbor preference by sample matrices for NEP method comparison

# load packages
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
# data path either do asymmetric or symmetric simulated data set
files = list.files("./../../../../data/20250217_sym00_nbh2_1000dim_grid200_300iter_50swaps/", pattern = ".csv")
data_path = "./../../../../data/20250217_sym00_nbh2_1000dim_grid200_300iter_50swaps/"
# output
output_path = "./../../../../../SCNA_thesis/github/Comparison/20250218_results_sym/histoCAT_sigval_delaunay_4ct_self00.csv"

### Data preparation
# load simulated data files into list
data_list = list()
for (i in files){
  data_list[[i]] <- readr::read_csv(paste0(data_path,i))
}
#add image id name
data_list = Map(cbind, data_list, img_id = names(data_list))
data = do.call("rbind", data_list)
# create dataframe with random values to pretend to have expression data to generate a giotto object
df = cbind(rep(1,base::nrow(data)),rep(1,base::nrow(data)))
colnames(df) = c("M1", "M2")
metadata = as.data.frame(cbind(as.character(data$ct), as.character(data$img_id), rownames(data)))
colnames(metadata) = c("ct", "img_id","cell_ID")

# Create Spatial Experiment
spe = SpatialExperiment(assay=as.data.frame(t(df)), 
                        colData=metadata, 
                        spatialCoords = as.matrix(data[, !(names(data) %in% c("ct", "img_id"))])
)

### NEP analysis
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = nbh_type,
                         k = k_number,
                         coords = c("x", "y"))

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
                          colPairName = "delaunay_interaction_graph",
                          BPPARAM = BPPARAM
                        )

### Create standard matrix output for comparison
#replace values that are not significant with False
#out$p_lt[out$sig == FALSE] <- 0.5

data = as.data.frame(out) %>% select(from_label, to_label, sigval, group_by)
data <- data.frame(lapply(data, function(x) ifelse(is.na(x) | is.nan(x), 0, x)))
data = as.data.frame(data)
data$key = paste(data$from_label, data$to_label, sep = "_")
data = data %>% select(-c(from_label, to_label))
data = spread(data, key = key, value = sigval)

rownames(data) = sub(".csv", "", data$group_by)
data = data[,-1]

write.csv(data,file=output_path,row.names = TRUE)

### Heatmap
library(pheatmap)
pheatmap(data)


sessionInfo()

writeLines(capture.output({
  print(R.version.string)  
  sessioninfo::session_info()  
}), "session_info.txt")


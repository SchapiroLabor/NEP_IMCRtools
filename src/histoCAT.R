# IMCR tools for simulated data

library(cytomapper)
library(imcRtools)
library(SingleCellExperiment)
library(readr)
library(SpatialExperiment)
library(tidyverse)

### This script runs IMCRtools histoCAT or classic for extracting interaction by sample matrices
# load in data

files = list.files("./../../../../data/Sim_nbh15_asym01_1000_grid0.2_1kiter_025kswap/", pattern = ".csv")
data_path = "./../../../../data/Sim_nbh15_asym01_1000_grid0.2_1kiter_025kswap/"
data_list = list()

for (i in files){
  data_list[[i]] <- readr::read_csv(paste0(data_path,i))
}


#add image id name
data_list = Map(cbind, data_list, img_id = names(data_list))
# add other metadata information and join all datasets

data = do.call("rbind", data_list)

# create dataframe with random values to pretend to have expression data to generate a giotto object
df = cbind(rep(1,base::nrow(data)),rep(1,base::nrow(data)))
#df <- as.matrix(runif(nrow(data) * 2), nrow = nrow(data))
colnames(df) = c("M1", "M2")

metadata = as.data.frame(cbind(as.character(data$ct), as.character(data$img_id), rownames(data)))
colnames(metadata) = c("ct", "img_id","cell_ID")

spe = SpatialExperiment(assay=as.data.frame(t(df)), 
                        colData=metadata, 
                        spatialCoords = as.matrix(data[, !(names(data) %in% c("ct", "img_id"))])
)

# delaunay, knn or distance neighborhood definition
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = "delaunay",
                         coords = c("x", "y"))

assayNames(spe) = "expr"

#adaps colPairName name depending on neighborhood definition

spe <- aggregateNeighbors(spe,
                          colPairName = "delaunay_interaction_graph",
                          aggregate_by = "metadata",
                          count_by = "ct") 

out <- countInteractions(spe,
                             group_by = "img_id",
                             label = "ct",
                             method = "classic",
                             colPairName = "delaunay_interaction_graph")

out <- testInteractions(spe,
                        group_by = "img_id",
                        label = "ct",
                        method = "classic",
                        colPairName = "delaunay_interaction_graph")

## create standard matrix for comparison
data = as.data.frame(out) %>% select(from_label, to_label, p_lt, group_by)
data$key = paste(data$from_label, data$to_label, sep = "_")
data = data %>% select(-c(from_label, to_label))
data = spread(data, key = key, value = p_lt)

rownames(data) = data$group_by
data = data[,-1]

write.csv(data,"./../../../Comparison/results_4ct_asym_0.2grid_self/IMCR_classic_delaunay_4ct_cross01.csv")

library(pheatmap)
pheatmap(data)

sessionInfo()


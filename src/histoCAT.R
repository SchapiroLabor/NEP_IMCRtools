# IMCR tools for simulated data

library(cytomapper)
library(imcRtools)
library(SingleCellExperiment)
library(readr)
library(SpatialExperiment)
library(tidyverse)
library(ggplot2)
library(ggraph)

# load in data
getwd()

files = list.files("./../../../../data/", pattern = ".csv")
data_path = "./../../../../data/"
data_list = list()

for (i in files){
  print(i)
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


spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = "expansion",
                         threshold = 20,
                         coords = c("x", "y"))
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = "knn",
                         k = 20,
                         coords = c("x", "y")) 
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = "delaunay",
                         coords = c("x", "y"))

assayNames(spe) = "expr"

spe <- aggregateNeighbors(spe,
                          colPairName = "knn_interaction_graph",
                          aggregate_by = "expression",
                          assay_type = "expr"
)

#head(spe$mean_aggregatedExpression)

spe <- aggregateNeighbors(spe,
                          colPairName = "delaunay_interaction_graph",
                          aggregate_by = "metadata",
                          count_by = "ct") 
head(spe$aggregatedNeighbors)

out <- countInteractions(spe,
                         group_by = "img_id",
                         label = "ct",
                         method = "histocat",
                         patch_size = 2,
                         colPairName = "delaunay_interaction_graph")

out <- testInteractions(spe,
                        group_by = "img_id",
                        label = "ct",
                        method = "histocat",
                        colPairName = "delaunay_interaction_graph")


write.csv(out,file=paste0("./../../output/histocat_output.csv"),row.names = TRUE)


#ggplot(as.data.frame(out), aes(from_label, to_label, fill = sigval)) +
#  facet_wrap(~ group_by, nrow = 6) +
#  geom_tile()


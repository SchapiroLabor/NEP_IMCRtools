# IMCR tools for simulated data

library(cytomapper)
library(imcRtools)
library(SingleCellExperiment)
library(readr)
library(SpatialExperiment)
library(tidyverse)


# load in data
getwd()

files = list.files("./../../../../data/Sim_100_asym_01/", pattern = ".csv")
data_path = "./../../../../data/Sim_100_asym_01/"
data_list = list()

for (i in files){
  data_list[[i]] <- readr::read_csv(paste0(data_path,i))
}

#sessionInfo()

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


#spe <- buildSpatialGraph(spe, img_id = "img_id",
#                         type = "expansion",
#                         threshold = 20,
#                         coords = c("x", "y"))
#spe <- buildSpatialGraph(spe, img_id = "img_id",
#                         type = "knn",
#                         k = 5,
#                         coords = c("x", "y")) 
spe <- buildSpatialGraph(spe, img_id = "img_id",
                         type = "delaunay",
                         coords = c("x", "y"))

assayNames(spe) = "expr"

#spe <- aggregateNeighbors(spe,
#                          colPairName = "knn_interaction_graph",
#                          aggregate_by = "expression",
#                          assay_type = "expr"
#)

#head(spe$mean_aggregatedExpression)

spe <- aggregateNeighbors(spe,
                          colPairName = "delaunay_interaction_graph",
                          aggregate_by = "metadata",
                          count_by = "ct") 

#spe <- aggregateNeighbors(spe,
 #                         colPairName = "knn_interaction_graph",
 #                         aggregate_by = "metadata",
 #                         count_by = "ct") 
#head(spe$aggregatedNeighbors)

#out_knn <- countInteractions(spe,
#                         group_by = "img_id",
#                         label = "ct",
#                         method = "histocat",
#                         #patch_size = 2,
#                         colPairName = "knn_interaction_graph")

out <- countInteractions(spe,
                             group_by = "img_id",
                             label = "ct",
                             method = "histocat",
                             #patch_size = 2,
                             colPairName = "delaunay_interaction_graph")

out <- testInteractions(spe,
                        group_by = "img_id",
                        label = "ct",
                        method = "histocat",
                        colPairName = "delaunay_interaction_graph")

#write.csv(out,file=paste0("./../../output/histocat_output_sim_100_01pref.csv"),row.names = TRUE)


## create standard matrix for comparison
data = as.data.frame(out) %>% select(from_label, to_label, p_lt, group_by)

data$key = paste(data$from_label, data$to_label, sep = "_")
data = data %>% select(-c(from_label, to_label))

data = spread(data, key = key, value = p_lt)


rownames(data) = data$group_by
data = data[,-1]

write.csv(data,"./../../../Comparison/results_4ct_cross01/IMCR_histoCAT_delaunay_p_lt_4ct_cross01.csv")


#ggplot(as.data.frame(out), aes(from_label, to_label, fill = sigval)) +
#  facet_wrap(~ group_by, nrow = 6) +
#  geom_tile()


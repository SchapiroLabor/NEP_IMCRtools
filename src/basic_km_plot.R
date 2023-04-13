### deal with output from count interaction without statistical testing
library(cytomapper)
library(imcRtools)
library(SingleCellExperiment)
library(readr)
library(SpatialExperiment)
library(tidyverse)
library(pheatmap)

out = out_del
out = out_knn
out = out

try = as.data.frame(out) %>% select(from_label, to_label, p_lt, group_by)
try = as.data.frame(out) %>% select(from_label, to_label, ct, group_by)
table(try$group_by)

### From here k-means, but rather try hierarchical clustering
try$key = paste(try$from_label, try$to_label, sep = "_")
try = try %>% select(-c(from_label, to_label))

try = spread(try, key = key, value = ct)
try = spread(try, key = key, value = p_lt)


library(ggplot2)
library(pheatmap)

rownames(try) = try$group_by
pheatmap(try[,-1], treeheight_row = 1, treeheight_col = 1)
pheatmap(try[,-1], kmeans_k = 2, treeheight_row = 0, 
         clustering_distance_cols = "euclidean",
         treeheight_col = 0)

pheatmap(try[,-1], clustering_distance_rows="euclidean", clustering_method="median", 
         kmeans_k = 2, show_rownames=TRUE)

km <- kmeans(try[,-1], centers=2)
cluster_assignments <- km$cluster
result = cbind(try[,1],cluster_assignments)

out = as.data.frame(out)

ggplot(out, aes(group_by, ct)) +
  boxplot()


### agglomerative clustering

set.seed(123)
# take try dataset from above and scale it
# normalize/scale
seeds_df_sc <- as.data.frame(scale(try[,-1]))
summary(seeds_df_sc)

# distance matrix
# distances between rows of dataset
dist_mat <- dist(seeds_df_sc, method = 'euclidean')

hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)












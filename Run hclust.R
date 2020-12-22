library(scater)
library(SC3)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(DuoClustering2018)
library(ggplot2)
library(mclust)
library(ClusterR)
library(cowplot)
library(reshape2)
library(cluster)
library(fpc)
library(clusterSim)
library(dbscan)
# read objects from PC:
scale.data.dist <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\scale.data.dist")
pca.data.dist <-readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data.dist")
df.distance.scale.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.scale.data")
df.distance.pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.pca.data")
pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data")
#

pca.data_20 <- pca.data[1:20,]

pca.data_20.dist <- dist(t(pca.data_20))
hclust_avg <- hclust(pca.data_20.dist, method = 'average')
cut_avg <- cutree(hclust_avg, k =8)
table(cut_avg)




hclust_avg <- hclust(pca.data.dist, method = 'average')
cut_avg <- cutree(hclust_avg, k =8)
table(cut_avg)



Run_hclust<- function (k = 8)
{
  hclust_avg <- hclust(pca.data.dist, method = 'average')
  #  cut the dendrogram in order to create the desired number of clusters.
  
  cut_avg <- cutree(hclust_avg, k =k)
  table(cut_avg)
  
  
  hclust_out <- data.frame(row.names = names(cut_avg), 
                           "labels" = gsub("^\\d+|\\d+$", "", names(cut_avg)),
                           "infered.label" = cut_avg)
  hclust_out$true.label <- unclass(as.factor(hclust_out$labels))
  return(hclust_out)
}


Calculate_IVM <- function(cluster_out)
{
  # Calculate Adjusted Rand Index
  ARI <- adjustedRandIndex(cluster_out$true.label, cluster_out$infered.label)
  
  # Reshape the distance matrices, add infered labels and true labels.
  df.distance.pca.data$Var1.infered.label = cluster_out$infered.label[match(df.distance.pca.data$Var1, row.names(cluster_out))]
  df.distance.pca.data$Var1.true.label = cluster_out$true.label[match(df.distance.pca.data$Var1, row.names(cluster_out))]
  df.distance.pca.data$Var2.infered.label = cluster_out$infered.label[match(df.distance.pca.data$Var2, row.names(cluster_out))]
  df.distance.pca.data$Var2.true.label = cluster_out$true.label[match(df.distance.pca.data$Var2, row.names(cluster_out))]
  
  df.distance.scale.data$Var1.infered.label = cluster_out$infered.label[match(df.distance.scale.data$Var1, row.names(cluster_out))]
  df.distance.scale.data$Var1.true.label = cluster_out$true.label[match(df.distance.scale.data$Var1, row.names(cluster_out))]
  df.distance.scale.data$Var2.infered.label = cluster_out$infered.label[match(df.distance.scale.data$Var2, row.names(cluster_out))]
  df.distance.scale.data$Var2.true.label = cluster_out$true.label[match(df.distance.scale.data$Var2, row.names(cluster_out))]
  
  # Calculate Dunn Index based on pca distance data
  # Separation is defined by maxium diameter of a cluster.
  # Cohesion is defined by minuim neighbour distance of two cluster.
  separation <- max(df.distance.pca.data[df.distance.pca.data$Var1.infered.label ==df.distance.pca.data$Var2.infered.label,]$value)
  cohesion <- min(df.distance.pca.data[df.distance.pca.data$Var1.infered.label !=df.distance.pca.data$Var2.infered.label,]$value)
  
  Dunn.Index <- cohesion/separation
  silhouette.index <- mean (silhouette(cluster_out$infered.label, pca.data.dist)[,3])
  CH <- cluster.stats(d = pca.data.dist,cluster_out$infered.label)$ch
  DB <- index.DB(x = t(pca.data), cl = cluster_out$infered.label, d=pca.data.dist, centrotypes="centroids", p=2, q=2)$DB
  
  out <- list("ARI" = ARI, 
              "Dunn.Index" = Dunn.Index, 
              "Infered.Cluster.Number" = max(cluster_out$infered.label), 
              "Mean.silhouette.index" = silhouette.index, 
              "Calinski.index" = CH, 
              "Davies.Bouldin.Index" = DB)
  
  # output infered label as a vector
  return(out)
}

IVM_baseline <- Calculate_IVM(Run_hclust(8))

out_put <- data.frame(IVM_baseline)

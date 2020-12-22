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


# manually do PCA
sce_zheng <- sce_full_Zhengmix8eq()

seurat.zheng <- CreateSeuratObject(counts = counts(sce_zheng))
seurat.zheng.lognorm <- NormalizeData(seurat.zheng)
seurat.zheng.feature <- FindVariableFeatures(object = seurat.zheng.lognorm)
all.genes <- rownames(seurat.zheng.lognorm)
seurat.zheng.scale <- ScaleData(seurat.zheng.feature, features = all.genes)

seurat.zheng.pca <- RunPCA(seurat.zheng.scale, features = VariableFeatures(object = seurat.zheng.scale))

#1. Distance matrix based on scaled data
scale.data <- GetAssayData(seurat.zheng.pca, slot = "scale.data")
euclidean.distance <- function(x, y) 
  return (sqrt(sum((x-y)^2)))
distance.scale.data <- CustomDistance(scale.data, euclidean.distance)
df.distance.scale.data<- melt(as.matrix(distance.scale.data))


#2. Distance matrix based on PCA data
pca.data <- Embeddings(seurat.zheng.pca,reduction = "pca")
pca.data <- t(pca.data)
euclidean.distance <- function(x, y) 
  return (sqrt(sum((x-y)^2)))
distance.pca.data <- CustomDistance(pca.data, euclidean.distance)
df.distance.pca.data <- melt(as.matrix(distance.pca.data))


# Calculate distance object use dist () function

scale.data.dist <- dist(t(GetAssayData(seurat.zheng.pca, slot = "scale.data")))
pca.data.dist <-dist(Embeddings(seurat.zheng.pca,reduction = "pca"))































scale.data.dist <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\scale.data.dist")
pca.data.dist <-readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data.dist")
df.distance.scale.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.scale.data")
df.distance.pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.pca.data")
pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data")

i = 40
dbscan::dbscan(t(pca.data), eps = i, minPts = 5)

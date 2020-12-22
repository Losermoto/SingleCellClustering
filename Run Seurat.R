# This is the pipeline for Seurat clustering
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
# Prepare for clustering
sce_zheng <- sce_full_Zhengmix8eq()
seurat_obj <- CreateSeuratObject(counts = counts(sce_zheng))
seurat_obj_lognorm <- NormalizeData(seurat_obj)
seurat_obj_featuresel <- FindVariableFeatures(seurat_obj_lognorm)
all.genes <- rownames(seurat_obj_featuresel)
seurat_obj_scale <- ScaleData(seurat_obj_featuresel, features = all.genes)
seurat_obj_pca <- RunPCA(seurat_obj_scale, features = VariableFeatures(object = seurat_obj_scale))

seurat_obj_jackstraw <- JackStraw(seurat_obj_pca, num.replicate = 100)
seurat_obj_score <- ScoreJackStraw(seurat_obj_jackstraw, dims = 1:20)
seurat_obj_neighb <- FindNeighbors(seurat_obj_score)
seurat_obj_cluster <- FindClusters(seurat_obj_neighb)


# Calculate distance matrix

#1. Distance matrix based on scaled data
scale.data <- GetAssayData(seurat_obj_cluster, slot = "scale.data")
euclidean.distance <- function(x, y) 
  return (sqrt(sum((x-y)^2)))
distance.scale.data <- CustomDistance(scale.data, euclidean.distance)
df.distance.scale.data<- melt(as.matrix(distance.scale.data))

#2. Distance matrix based on PCA data
pca.data <- Embeddings(seurat_obj_cluster,reduction = "pca")
pca.data <- t(pca.data)
euclidean.distance <- function(x, y) 
  return (sqrt(sum((x-y)^2)))
distance.pca.data <- CustomDistance(pca.data, euclidean.distance)
df.distance.pca.data <- melt(as.matrix(distance.pca.data))


# Calculate distance object use dist () function

scale.data.dist <- dist(t(GetAssayData(seurat_obj_cluster, slot = "scale.data")))
pca.data.dist <-dist(Embeddings(seurat_obj_cluster,reduction = "pca"))

# Save two distance matrices and two dist objects
# saveRDS(pca.data, "D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data")
# saveRDS(scale.data.dist, "D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\scale.data.dist")
# saveRDS(pca.data.dist, "D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data.dist")

# saveRDS(df.distance.scale.data, "D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.scale.data")
# saveRDS(df.distance.pca.data, "D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.pca.data")

# scale.data.dist <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\scale.data.dist")
# pca.data.dist <-readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data.dist")
# 
# df.distance.scale.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.scale.data")
# df.distance.pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.pca.data")

# Function for clustering 
Run_seurat <- function(resolution){
  seurat_obj_cluster <- FindClusters(seurat_obj_neighb, resolution = resolution)
  seurat_out <- data.frame(Idents(seurat_obj_cluster))
  seurat_out$labels = gsub("^\\d+|\\d+$", "", row.names(seurat_out))
  seurat_out$true.label <- unclass(as.factor(seurat_out$labels))
  seurat_out$infered.label <- as.numeric(seurat_out$Idents.seurat_obj_cluster.)
  return(seurat_out)
}

# Function for calculaing IVM and ARI
# Input cluster_out
Calculate_IVM <- function(cluster_out)
{
  # Calculate Adjusted Rand Index
  ARI <- adjustedRandIndex(cluster_out$Idents.seurat_obj_cluster., cluster_out$labels)
  
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

  out <- list("ARI" = ARI, "Dunn.Index" = Dunn.Index, "Infered.Cluster.Number" = max(cluster_out$infered.label), 
              "Mean.silhouette.index" = silhouette.index, "Calinski.index" = CH, "Davies.Bouldin.Index" = DB)
  
  # output infered label as a vector
  return(out)

}
# Now run the function in the range of resolution (0.4-1.4)

# input the number of cycles


ARI.vector <- vector(mode = "numeric", length = 40)
Dunn.Index.vector <- vector(mode = "numeric", length = 40)
Infered.Cluster.Number.vector <- vector(mode = "numeric", length = 40)
silhouette.index.vector <- vector(mode = "numeric", length = 40)
Calinski.index.vector <- vector(mode = "numeric", length = 40)
Davies.Bouldin.Index.vector <-  vector(mode = "numeric", length = 40)
for (i in 1:40) {
  temp.out <- Run_seurat(0.025*i+.1)
  temp.out.2 <- Calculate_IVM(temp.out)
  ARI.vector[i] <- temp.out.2[[1]]
  Dunn.Index.vector[i] <- temp.out.2[[2]]
  Infered.Cluster.Number.vector[i]<- temp.out.2[[3]]
  silhouette.index.vector[i] <- temp.out.2[[4]]
  Calinski.index.vector [i] <- temp.out.2[[5]]
  Davies.Bouldin.Index.vector [i] <- temp.out.2[[6]]
}
x <- 1:40*0.025+0.1

export.df <- data.frame("resolution" = x,ARI.vector, Dunn.Index.vector, Infered.Cluster.Number.vector, 
                        silhouette.index.vector,Calinski.index.vector,Davies.Bouldin.Index.vector)



cor.test(export.df$ARI.vector, export.df$Dunn.Index.vector)
cor.test(export.df$ARI.vector, export.df$silhouette.index.vector)
cor.test(export.df$ARI.vector, export.df$Calinski.index.vector)
cor.test(export.df$ARI.vector, export.df$Davies.Bouldin.Index.vector)



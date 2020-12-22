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
# Prepare for SC3 clustering
sce_zheng_sc3 <- sce_full_Zhengmix8eq()
rowData(sce_zheng_sc3)$feature_symbol <- rownames(counts(sce_zheng_sc3))
sc3_object <- sc3_prepare(sce_zheng_sc3)
sc3_dist <- sc3_calc_dists(sc3_object)
sc3_transf <- sc3_calc_transfs(sc3_dist)
# sc3_cluster <- sc3_kmeans(sc3_transf, ks = 8)
# sc3_consensus <-sc3_calc_consens(sc3_cluster)

# saveRDS(sc3_transf, "D:\\New Term\\scRNA_clustering\\RNAClustering\\SC3\\sc3_transf")


sc3_transf <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\SC3\\sc3_transf")


scale.data.dist <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\scale.data.dist")
pca.data.dist <-readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data.dist")

df.distance.scale.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.scale.data")
df.distance.pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\df.distance.pca.data")

pca.data <- readRDS("D:\\New Term\\scRNA_clustering\\RNAClustering\\Seurat\\pca.data")


# Function for SC3:
t.start <- proc.time()

Run_SC3 <- function(k){
  sc3_cluster <- sc3_kmeans(sc3_transf, ks = k)
  sc3_consensus <-sc3_calc_consens(sc3_cluster)
  col_data <- colData(sc3_consensus)
  SC3_out <- data.frame(row.names = row.names(col_data), 
                        "labels" = col_data$phenoid, 
                        "infered.label" = as.numeric(col_data[,grep("sc3_", colnames(col_data))]))
  SC3_out$true.label <- unclass(as.factor(SC3_out$labels))
  return(SC3_out)
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
ARI.vector <- vector(mode = "numeric", length = 16)
Dunn.Index.vector <- vector(mode = "numeric", length = 16)
Infered.Cluster.Number.vector <- vector(mode = "numeric", length = 16)
silhouette.index.vector <- vector(mode = "numeric", length = 16)
Calinski.index.vector <- vector(mode = "numeric", length = 16)
Davies.Bouldin.Index.vector <-  vector(mode = "numeric", length = 16)


for (i in 4:19){
  temp.out<- Run_SC3( k = i)
  
  temp.out.2 <- Calculate_IVM(temp.out)
  ARI.vector[i] <- temp.out.2[[1]]
  Dunn.Index.vector[i] <- temp.out.2[[2]]
  Infered.Cluster.Number.vector[i]<- temp.out.2[[3]]
  silhouette.index.vector[i] <- temp.out.2[[4]]
  Calinski.index.vector [i] <- temp.out.2[[5]]
  Davies.Bouldin.Index.vector [i] <- temp.out.2[[6]]
}

x <-1:19

export.df <- data.frame("Number of cluster" = x,ARI.vector, Dunn.Index.vector, Infered.Cluster.Number.vector, 
                        silhouette.index.vector,Calinski.index.vector,Davies.Bouldin.Index.vector)

export.df <- export.df[-c(1,2,3),]
write.csv(export.df, "D:\\New Term\\scRNA_clustering\\RNAClustering\\SC3\\IVMs.csv")
# export.df <- read.csv("D:\\New Term\\scRNA_clustering\\RNAClustering\\SC3\\IVMs.csv")



time.elasped <- proc.time()-t.start
time.elasped



cor.test(export.df$ARI.vector, export.df$Dunn.Index.vector)
cor.test(export.df$ARI.vector, export.df$silhouette.index.vector)
cor.test(export.df$ARI.vector, export.df$Calinski.index.vector)
cor.test(export.df$ARI.vector, export.df$Davies.Bouldin.Index.vector)

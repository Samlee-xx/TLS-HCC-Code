
###### HD bile duct similarity change along axis ######
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggplot2)

bd.data <- readRDS('...\\Merge Normal HCC bile duct.rds')
localdir <- '...\\outs'
hd.data <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

bile.duct <- read.csv('...\\Bile Duct Cluster.csv')
colnames(hd.data)

hd.bile <- hd.data[,bile.duct$Barcode]
dim(hd.bile)
dim(hd.data)

barcode.xy <- read.csv('...\\Sample5 Barcode centroids.csv')
rownames(barcode.xy) = barcode.xy$barcode
barcode.xy <- barcode.xy[bile.duct$Barcode,]

a=ggplot(barcode.xy)+
  geom_point(aes(x,-y),size=0.5)+
  theme_classic()

pdf(file = '...\\Bileduct distribution.pdf',
    height = 5,width = 5)
print(a)
dev.off()

expr_matrix <- hd.bile@assays$Spatial.008um@layers$counts %>% as.matrix()
colnames(expr_matrix) = colnames(hd.bile)
rownames(expr_matrix) = rownames(hd.bile)

dim(expr_matrix)
rownames(expr_matrix)
sum(!colnames(hd.bile) == barcode.xy$barcode)
barcode_data <- data.frame(Barcode = colnames(expr_matrix), y = barcode.xy$y)

used.barcode <- filter(barcode.xy,x < 7500)
common.genes <- intersect(rownames(bd.data), rownames(expr_matrix))
hcc.avg <- rowMeans(bd.data@assays$RNA@data[common.genes, bd.data$group == "hcc"])
health.avg <- rowMeans(bd.data@assays$RNA@data[common.genes, bd.data$group == "health"])
hd.expr <- expr_matrix[common.genes, ] 

cosine_similarity <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}

sim.hcc <- numeric(ncol(hd.expr))
sim.health <- numeric(ncol(hd.expr))

for (i in 1:ncol(hd.expr)) {
  spot.vec <- hd.expr[, i]
  sim.hcc[i] <- cosine_similarity(spot.vec, hcc.avg)
  sim.health[i] <- cosine_similarity(spot.vec, health.avg)
}

similarity_df <- data.frame(
  barcode = colnames(hd.expr),
  x = barcode.xy$x,
  sim_hcc = sim.hcc,
  sim_health = sim.health
)

rownames(similarity_df) <- similarity_df$barcode
similarity_df <- similarity_df[used.barcode$barcode,]
similarity_df$sim_hcc_scale <- scale(similarity_df$sim_hcc)
similarity_df$sim_health_scale <- scale(similarity_df$sim_health)
similarity_df$diff <- similarity_df$sim_hcc_scale - similarity_df$sim_health_scale

a = ggplot(similarity_df, aes(x = x, y = diff)) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", size = 1.2) +
  theme_classic() +
  labs(
    title = "Trend of Tumor Similarity Along Y Axis",
    y = "HCC similarity - Health similarity",
    x = "Y axis (proximity to tumor)"
  )

pdf(file = '...\\Similarity2.pdf',
    height = 3,width = 3)
print(a)
dev.off()



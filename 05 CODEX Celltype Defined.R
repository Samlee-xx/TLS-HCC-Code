
setwd('...\\09 Celltype Defined')
library(Seurat)
library(reshape2)
library(dplyr)
tls = readRDS('Seurat TLS.rds')
celltype.marker = c('APOC3','COX2','Pan.Cytokeratin', #tumor
                    'MECA79','CD31', 'Podoplanin',# Endothelial
                    'aSMA', # Fibroblast
                    'CD45','CD45RO', # Immune marker
                    'CD14','CD66b','CD68','CD11c','MPO', # Myeloid
                    'CD20', 'CD21','CD79a', # B cell
                    'CD3e','CD56','CD4','CD8'  # NKT 
)

a=DoHeatmap(subset(tls, downsample=3000),features = celltype.marker,angle = 90, size = 4,
            disp.min = -2, disp.max = 2)

pdf(file = '01 DoHeatmap Celltype marker.pdf',height = 8,width = 20)
print(a)
dev.off()

feature.filtered = read.csv('...\\Merge_total.csv')
colnames(feature.filtered)
feature.filtered = feature.filtered[,c(1,2,53,54,55,56,57)]
colnames(feature.filtered)[1] = 'cellid'
rownames(feature.filtered) = feature.filtered$cellid
feature.filtered = feature.filtered[colnames(tls),]
saveRDS(feature.filtered,'Filter size Metadta.rds')

cell.metadata.total = readRDS('Filter size Metadta.rds')
dim(cell.metadata.total)
sum(!colnames(tls) == cell.metadata.total$cellid)
tls$fov = cell.metadata.total$fov

df.fov = table(tls$seurat_clusters,tls$fov) %>% as.data.frame()
df.fov = df.fov %>% group_by(Var1) %>% top_n(10,Freq)

for (i in c(0:50)){
  check.cluster = subset(tls,seurat_clusters %in% i)
  check.cluster$celltype = paste0('Cluster',i)
  
  check.metadata = data.frame(cellid = colnames(check.cluster),celltype = check.cluster$celltype)
  check.metadata2 = check.metadata %>% left_join(cell.metadata.total)
  
  output_filename <- paste0('Check_cluster', i, '.csv') 
  output_path <- file.path('...\\check clusters', output_filename) 
  write.csv(check.metadata2, output_path, row.names = FALSE)
}

df.fov = filter(df.fov,Freq >=50)
write.csv(df.fov,'fov information check.csv')

library(stringr)
df.fov2 = str_split(df.fov$Var2,'_',simplify = TRUE)
df.fov2 = df.fov2 %>% as.data.frame()
df.fov$slide = paste0('Sample_',df.fov2$V1)
table(df.fov$slide)
write.csv(df.fov,'fov information check.csv')

immune.cluster = c(0,1,2,4,5,7,15,20,22,27,32,41,44,48)
stroma.cluster = c(6,11,16)
tumor.cluster = c(8,9,10,13,14,17,18,19,21,23,24,25,26,28,30,31,33,36,37,38,39,
                  40,42,43,46,47,49)
artifact.cluster = c(12,34,45,50)
mix.cluster = c(3,29,35)

immune.cluster = subset(tls,seurat_clusters %in% immune.cluster)
stroma.cluster = subset(tls,seurat_clusters %in% stroma.cluster)
tumor.cluster = subset(tls,seurat_clusters %in% tumor.cluster)
artifact.cluster = subset(tls,seurat_clusters %in% artifact.cluster)

immune.cluster$celltype = 'Immune'
stroma.cluster$celltype = 'Stroma'
tumor.cluster$celltype = 'Tumor'
artifact.cluster$celltype = 'Artifact'

mix.cluster = subset(tls,seurat_clusters %in% mix.cluster)
VariableFeatures(mix.cluster) <- celltype.marker
mix.cluster = ScaleData(mix.cluster)
mix.cluster <- RunPCA(object = mix.cluster, npcs = 20, verbose = TRUE)
mix.cluster <- RunUMAP(object = mix.cluster, dims = 1:20, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
mix.cluster <- FindNeighbors(object = mix.cluster, dims = 1:20, verbose = TRUE)
mix.cluster <- FindClusters(object = mix.cluster, verbose = TRUE, resolution = 0.3)
DimPlot(mix.cluster)
table(mix.cluster$seurat_clusters)
saveRDS(mix.cluster,'Seurat Mix.rds')
mix.cluster = readRDS('Seurat Mix.rds')

df.fov = table(mix.cluster$seurat_clusters,mix.cluster$fov) %>% as.data.frame()
df.fov = df.fov %>% group_by(Var1) %>% top_n(10,Freq)

for (i in c(0:13)){
  check.cluster = subset(mix.cluster,seurat_clusters %in% i)
  check.cluster$celltype = paste0('Cluster',i)
  
  check.metadata = data.frame(cellid = colnames(check.cluster),celltype = check.cluster$celltype)
  check.metadata2 = check.metadata %>% left_join(cell.metadata.total)
  
  output_filename <- paste0('Check_cluster', i, '.csv') 
  output_path <- file.path('...\\check mix clusters', output_filename) 
  write.csv(check.metadata2, output_path, row.names = FALSE)
}

df.fov = filter(df.fov,Freq >=100)
library(stringr)
df.fov2 = str_split(df.fov$Var2,'_',simplify = TRUE)
df.fov2 = df.fov2 %>% as.data.frame()
df.fov$slide = paste0('Sample_',df.fov2$V1)
table(df.fov$slide)
write.csv(df.fov,'fov mix information check.csv')

a=DoHeatmap(subset(mix.cluster, downsample=3000),features = celltype.marker,angle = 90, size = 4,
            disp.min = -2, disp.max = 2)

pdf(file = '02 DoHeatmap mix cluster Celltype marker.pdf',height = 8,width = 20)
print(a)
dev.off()

mix.immune = c(0,3,4,5)
mix.tumor = c(6,8,9,10,11,12)
mix.stroma = c(1,13)
mix.artifact = c(2)
mix.undefined = c(7)

mix.immune = subset(mix.cluster,seurat_clusters %in% mix.immune)
mix.tumor = subset(mix.cluster,seurat_clusters %in% mix.tumor)
mix.stroma = subset(mix.cluster,seurat_clusters %in% mix.stroma)
mix.artifact = subset(mix.cluster,seurat_clusters %in% mix.artifact)
mix.undefined = subset(mix.cluster,seurat_clusters %in% mix.undefined)

mix.immune$celltype = 'Immune'
mix.tumor$celltype = 'Tumor'
mix.stroma$celltype = 'Stroma'
mix.artifact$celltype = 'Artifact'
mix.undefined$celltype = 'Undefined'

celltype.metadata = data.frame(cellid = c(colnames(immune.cluster),
                                          colnames(stroma.cluster),
                                          colnames(tumor.cluster),
                                          colnames(artifact.cluster),
                                          colnames(mix.immune),
                                          colnames(mix.tumor),
                                          colnames(mix.stroma),
                                          colnames(mix.artifact),
                                          colnames(mix.undefined)),
                               celltype = c(immune.cluster$celltype,
                                            stroma.cluster$celltype,
                                            tumor.cluster$celltype,
                                            artifact.cluster$celltype,
                                            mix.immune$celltype,
                                            mix.tumor$celltype,
                                            mix.stroma$celltype,
                                            mix.artifact$celltype,
                                            mix.undefined$celltype))

cell.metadata.total = readRDS('Filter size Metadta.rds')
celltype.metadata = celltype.metadata %>% left_join(cell.metadata.total)

celltype.metadata$slide = ifelse(celltype.metadata$slide == 'sample_136759','Sample_136759',celltype.metadata$slide)

write.csv(celltype.metadata,'CT defined metadata.csv')
saveRDS(celltype.metadata,'CT defined metadata.rds')

rownames(celltype.metadata) = celltype.metadata$cellid
celltype.metadata = celltype.metadata[colnames(tls),]
tls$celltype = celltype.metadata$celltype
DimPlot(tls,split.by  = 'celltype')











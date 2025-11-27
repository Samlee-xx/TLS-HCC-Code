
###### Task 2 ######

library(dplyr)
library(reshape2)
library(stringr)
library(sva)
library(Seurat)
setwd('...\\Analysis')
data = read.csv('...\\cell_table_arcsinh_transformed.csv')
data.sub = filter(data, cell_size <=300)
data.sub = filter(data.sub, cell_size >=30)
rownames(data.sub) = paste0(data.sub$fov,"_",data.sub$label)

count.data = data.sub[,c(2:15,17:34)]
rownames(count.data) = paste0(data.sub$fov,"_",data.sub$label)
count.data = t(count.data)

data.obj = CreateSeuratObject(counts =  count.data)
colnames(data.obj)
rownames(data.obj)
celltype.marker = c('EpCAM','Pan.Cytokeratin','Hepar1', #tumor
                    'CD31','Podoplanin',# Endothelial
                    'SMA',
                    'CD45','CD3e','CD20','CD4','CD8','CD68','MPO','S100A8_9' # Immune marker
)

VariableFeatures(data.obj) <- rownames(data.obj)
data.obj = ScaleData(data.obj)
data.obj@assays$RNA@scale.data[1:15,1:5]
data.obj <- RunPCA(object = data.obj, npcs = 10, verbose = TRUE)
data.obj <- RunUMAP(object = data.obj, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
data.obj <- FindNeighbors(object = data.obj, dims = 1:10, verbose = TRUE)
data.obj <- FindClusters(object = data.obj, verbose = TRUE, resolution = 0.8)
DimPlot(data.obj)

a = DoHeatmap(subset(data.obj, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '01 DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.cluster = subset(data.obj,seurat_clusters %in% c(3,4,6,7,9,16,19,20,21))
epi.cluster = subset(data.obj,seurat_clusters %in% c(0,1,2,5,8,12,13,14,15))
stroma.cluster = subset(data.obj,seurat_clusters %in% c(11))
unassign.cluster = subset(data.obj,seurat_clusters %in% c(10,17,18,22,23))

rownames(immune.cluster)
immune.marker = c('CD14','CD163','CD68','CD20','CD79a','CD3e','CD4','CD8','S100A8-9','MPO')
VariableFeatures(immune.cluster) <- immune.marker
immune.cluster = ScaleData(immune.cluster)
immune.cluster <- RunPCA(object = immune.cluster, npcs = 5, verbose = TRUE)
immune.cluster <- RunUMAP(object = immune.cluster, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.cluster <- FindNeighbors(object = immune.cluster, dims = 1:5, verbose = TRUE)
immune.cluster <- FindClusters(object = immune.cluster, verbose = TRUE, resolution = 1)
DimPlot(immune.cluster)
rownames(immune.cluster)

lymphoid = subset(immune.cluster,seurat_clusters %in% c(11,13,3,5,14,18,23))
myeloid = subset(immune.cluster,seurat_clusters %in% c(7,9,19,12,20,22))
immune.unassign = subset(immune.cluster,seurat_clusters %in% c(21))

immune.secend = subset(immune.cluster,seurat_clusters %in% c(0,1,2,4,6,8,10,15,16,17))
a = DoHeatmap(subset(immune.cluster, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '01 DoHeatmap IMMUNE marker purity.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.marker = c('CD14','CD163','CD68','CD20','CD79a','CD3e','CD4','CD8','S100A8-9','MPO')
VariableFeatures(immune.secend) <- immune.marker
immune.secend = ScaleData(immune.secend)
immune.secend <- RunPCA(object = immune.secend, npcs = 5, verbose = TRUE)
immune.secend <- RunUMAP(object = immune.secend, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.secend <- FindNeighbors(object = immune.secend, dims = 1:5, verbose = TRUE)
immune.secend <- FindClusters(object = immune.secend, verbose = TRUE, resolution = 1)
DimPlot(immune.secend)
rownames(immune.secend)

a = DoHeatmap(subset(immune.secend, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '02 DoHeatmap immune secend.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.second.lymphoid = subset(immune.secend,seurat_clusters %in% c(0,2,4,5,12,14,15,18,21,22,23))
immune.second.myeloid = subset(immune.secend,seurat_clusters %in% c(1,3,6,7,8,9,10,17,19,20))
immune.second.tumor =  subset(immune.secend,seurat_clusters %in% c(13))
#11 16 removed

second.cluster.file = data.sub[colnames(immune.second.check),]
second.cluster.file$cluster =paste0("Cluster_",immune.second.check$seurat_clusters) 
write.csv(second.cluster.file,'...\\check1.csv')

# lymphoid
lymphoid.total = data.obj[,c(colnames(lymphoid),colnames(immune.second.lymphoid))]
rownames(lymphoid.total)
lymphoid.marker = c('CD3e','CD4','CD8','CD20','CD79a')
VariableFeatures(lymphoid.total) <- lymphoid.marker
lymphoid.total = ScaleData(lymphoid.total)
lymphoid.total <- RunPCA(object = lymphoid.total, npcs = 3, verbose = TRUE)
lymphoid.total <- RunUMAP(object = lymphoid.total, dims = 1:3, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
lymphoid.total <- FindNeighbors(object = lymphoid.total, dims = 1:3, verbose = TRUE)
lymphoid.total <- FindClusters(object = lymphoid.total, verbose = TRUE, resolution = 0.5)
DimPlot(lymphoid.total)
rownames(immune.secend)
a = DoHeatmap(subset(lymphoid.total, downsample=2000), features = lymphoid.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)

cd8t = subset(lymphoid.total,seurat_clusters %in% c(0,1,5,8))
cd4t = subset(lymphoid.total,seurat_clusters %in% c(3,7,13,14))
t = subset(lymphoid.total,seurat_clusters %in% c(12,17))
b = subset(lymphoid.total,seurat_clusters %in% c(2,4,6,9,11,18,20))
t.b = subset(lymphoid.total,seurat_clusters %in% c(10,21))

# myeloid
myeloid.total = data.obj[,c(colnames(myeloid),colnames(immune.second.myeloid))]
rownames(myeloid.total)
myeloid.marker = c('CD14','CD68','CD163','MPO','S100A8-9')
VariableFeatures(myeloid.total) <- myeloid.marker
myeloid.total = ScaleData(myeloid.total)
myeloid.total <- RunPCA(object = myeloid.total, npcs = 3, verbose = TRUE)
myeloid.total <- RunUMAP(object = myeloid.total, dims = 1:3, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
myeloid.total <- FindNeighbors(object = myeloid.total, dims = 1:3, verbose = TRUE)
myeloid.total <- FindClusters(object = myeloid.total, verbose = TRUE, resolution = 0.5)
DimPlot(myeloid.total)
rownames(myeloid.total)
a = DoHeatmap(subset(myeloid.total, downsample=2000), features = myeloid.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)

macrophage = subset(myeloid.total,seurat_clusters %in% c(0,2,4,6,9,10,12,13,14,15,16))
neutrophil = subset(myeloid.total,seurat_clusters %in% c(1,3,5,7,8,11))

# tumor 
epi.cluster = subset(data.obj,seurat_clusters %in% c(0,1,2,5,8,12,13,14,15))
rownames(epi.cluster)
tumor.marker = c('EpCAM','Hepar1','Pan.Cytokeratin')
VariableFeatures(epi.cluster) <- tumor.marker
epi.cluster = ScaleData(epi.cluster)
epi.cluster <- RunPCA(object = epi.cluster, npcs = 2, verbose = TRUE)
epi.cluster <- RunUMAP(object = epi.cluster, dims = 1:2, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
epi.cluster <- FindNeighbors(object = epi.cluster, dims = 1:2, verbose = TRUE)
epi.cluster <- FindClusters(object = epi.cluster, verbose = TRUE, resolution = 0.3)
a = DoHeatmap(subset(epi.cluster, downsample=2000), features = tumor.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)

bile.duct = subset(epi.cluster,seurat_clusters %in% 2)
tumor = subset(epi.cluster,seurat_clusters %in% c(0,1,3:43))

##total
bile.duct = subset(epi.cluster,seurat_clusters %in% 2)
tumor = subset(epi.cluster,seurat_clusters %in% c(0,1,3:43))

cd8t = subset(lymphoid.total,seurat_clusters %in% c(0,1,5,8))
cd4t = subset(lymphoid.total,seurat_clusters %in% c(3,7,13,14))
t = subset(lymphoid.total,seurat_clusters %in% c(12,17))
b = subset(lymphoid.total,seurat_clusters %in% c(2,4,6,9,11,18,20))
t.b = subset(lymphoid.total,seurat_clusters %in% c(10,21))
macrophage = subset(myeloid.total,seurat_clusters %in% c(0,2,4,6,9,10,12,13,14,15,16))
neutrophil = subset(myeloid.total,seurat_clusters %in% c(1,3,5,7,8,11))
stroma.cluster = subset(data.obj,seurat_clusters %in% c(11))

data.sub[colnames(stroma.cluster),'celltype'] = 'Endothelial'
data.sub[colnames(neutrophil),'celltype'] = 'Neutrophil'
data.sub[colnames(t.b),'celltype'] = 'T-B'
data.sub[colnames(b),'celltype'] = 'B_cell'
data.sub[colnames(t),'celltype'] = 'T_cell'
data.sub[colnames(cd4t),'celltype'] = 'CD4T'
data.sub[colnames(cd8t),'celltype'] = 'CD8T'
data.sub[colnames(tumor),'celltype'] = 'Tumor'
data.sub[colnames(bile.duct),'celltype'] = 'Bile_Duct'
data.sub[colnames(macrophage),'celltype'] = 'Macrophage'
data.sub[is.na(data.sub$celltype),'celltype'] = 'Unassign'
table(data.sub$celltype)
write.csv(data.sub,'...\\Add celltype.csv')

# advanvced 
#T-B need to divided; Bile Duct need to purity
t.b = subset(lymphoid.total,seurat_clusters %in% c(10,21))
lymphoid.marker = c('CD3e','CD4','CD8','CD20','CD79a')
VariableFeatures(t.b) <- lymphoid.marker
t.b = ScaleData(t.b)
t.b <- RunPCA(object = t.b, npcs = 3, verbose = TRUE)
t.b <- RunUMAP(object = t.b, dims = 1:3, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
t.b <- FindNeighbors(object = t.b, dims = 1:3, verbose = TRUE)
t.b <- FindClusters(object = t.b, verbose = TRUE, resolution = 0.5)
DimPlot(t.b)
DoHeatmap(subset(t.b, downsample=2000), features = lymphoid.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)

t.b.b = subset(t.b,seurat_clusters %in% c(2,6))
t.b.t = subset(t.b,seurat_clusters %in% c(4,5,7))
t.b.cd4t = subset(t.b,seurat_clusters %in% c(0))
t.b.cd8t = subset(t.b,seurat_clusters %in% c(1,3))

# Bile duct 
rownames(bile.duct)
tumor.marker = c('EpCAM','Hepar1','Pan.Cytokeratin')
VariableFeatures(bile.duct) <- tumor.marker
bile.duct = ScaleData(bile.duct)
bile.duct <- RunPCA(object = bile.duct, npcs = 2, verbose = TRUE)
bile.duct <- RunUMAP(object = bile.duct, dims = 1:2, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
bile.duct <- FindNeighbors(object = bile.duct, dims = 1:2, verbose = TRUE)
bile.duct <- FindClusters(object = bile.duct, verbose = TRUE, resolution = 0.3)
DimPlot(bile.duct)
table(bile.duct$seurat_clusters)
DoHeatmap(subset(bile.duct, downsample=2000), features = tumor.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)

bile.duct.file = data.sub[colnames(bile.duct),]
bile.duct.file$cluster =paste0("Cluster_",bile.duct$seurat_clusters) 
write.csv(bile.duct.file,'...\\checkBileDuct.csv')

bile.duct.bd =subset(bile.duct,seurat_clusters %in% c(2,4,7,11,13))
bile.duct.tumor =subset(bile.duct,seurat_clusters %in% c(0,1,3,5,6,8,9,10,12,14))

# total final
t.b.b = subset(t.b,seurat_clusters %in% c(2,6))
t.b.t = subset(t.b,seurat_clusters %in% c(4,5,7))
t.b.cd4t = subset(t.b,seurat_clusters %in% c(0))
t.b.cd8t = subset(t.b,seurat_clusters %in% c(1,3))

data.sub[colnames(stroma.cluster),'celltype'] = 'Endothelial'
data.sub[colnames(neutrophil),'celltype'] = 'Neutrophil'
#data.sub[colnames(t.b),'celltype'] = 'T-B'
data.sub[c(colnames(b),colnames(t.b.b)),'celltype'] = 'B_cell'
data.sub[c(colnames(t),colnames(t.b.t)),'celltype'] = 'T_cell'
data.sub[c(colnames(cd4t),colnames(t.b.cd4t)),'celltype'] = 'CD4T'
data.sub[c(colnames(cd8t),colnames(t.b.cd8t)),'celltype'] = 'CD8T'
data.sub[c(colnames(tumor),colnames(bile.duct.tumor)),'celltype'] = 'Tumor'
data.sub[colnames(bile.duct.bd),'celltype'] = 'Bile_Duct'
data.sub[colnames(macrophage),'celltype'] = 'Macrophage'
data.sub[is.na(data.sub$celltype),'celltype'] = 'Unassign'
write.csv(data.sub,'...\\Add celltype.csv')
table(data.sub$celltype)

###### Task 1 ######


setwd('...\\Analysis')
data = read.csv('...\\cell_table_arcsinh_transformed.csv')

data.sub = filter(data, cell_size <=300)
data.sub = filter(data.sub, cell_size >=30)
rownames(data.sub) = paste0(data.sub$fov,"_",data.sub$label)

colnames(data.sub)
count.data = data.sub[,c(2:15,17:34)]
rownames(count.data) = paste0(data.sub$fov,"_",data.sub$label)
count.data = t(count.data)

data.obj = CreateSeuratObject(counts =  count.data)
colnames(data.obj)
rownames(data.obj)
celltype.marker = c('EpCAM','Pan.Cytokeratin','Hepar1', #tumor
                    'CD31','Podoplanin',# Endothelial
                    'SMA',
                    'CD45','CD3e','CD20','CD4','CD8','CD68','MPO' # Immune marker
)

VariableFeatures(data.obj) <- rownames(data.obj)
data.obj = ScaleData(data.obj)
data.obj@assays$RNA@scale.data[1:15,1:5]
data.obj <- RunPCA(object = data.obj, npcs = 10, verbose = TRUE)
data.obj <- RunUMAP(object = data.obj, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
data.obj <- FindNeighbors(object = data.obj, dims = 1:10, verbose = TRUE)
data.obj <- FindClusters(object = data.obj, verbose = TRUE, resolution = 0.8)
DimPlot(data.obj)

a = DoHeatmap(subset(data.obj, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '01 DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

immune = subset(data.obj,seurat_clusters %in% c(0,1,3,7,10,11,15,16,17,18,20,22,23,26,27))
tumor = subset(data.obj,seurat_clusters %in% c(2,4,5,6,12,14,19,21,24,25,29))
mix = subset(data.obj,seurat_clusters %in% c(8,9,13))

# mix divided
VariableFeatures(mix) <- rownames(mix)
mix = ScaleData(mix)
mix@assays$RNA@scale.data[1:15,1:5]
mix <- RunPCA(object = mix, npcs = 10, verbose = TRUE)
mix <- RunUMAP(object = mix, dims = 1:10, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
mix <- FindNeighbors(object = mix, dims = 1:10, verbose = TRUE)
mix <- FindClusters(object = mix, verbose = TRUE, resolution = 0.8)
DimPlot(mix)

a = DoHeatmap(subset(mix, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '02 mix DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

mix.file = data.sub[colnames(mix),]
mix.file.part1 = subset(mix.file,cluster %in% c(0:8))
mix.file.part1$cluster =paste0("Cluster_",mix.file.part1$seurat_clusters) 

mix.fil.part1 = subset(mix,seurat_clusters %in% c(0:8))
mix.file.file.part1 = data.sub[colnames(mix.fil.part1),]
mix.file.file.part1$cluster =paste0("Cluster_",mix.fil.part1$seurat_clusters) 
write.csv(mix.file.file.part1,'Mix Part1.csv')

mix.immune = subset(mix,seurat_clusters %in% c(0,1,4,5,6))
mix.tumor = subset(mix,seurat_clusters %in% c(2,3,7,8,9,10,11,12,13,14,16))

# Immune and stroma Part
immune.part = data.obj[,c(colnames(immune,mix.immune))]
celltype.marker = c('CD31','Podoplanin',# Endothelial
                    'SMA',
                    'CD45','CD3e','CD20','CD79a','CD68','MPO' # Immune marker
)
VariableFeatures(immune.part) <- celltype.marker
immune.part = ScaleData(immune.part)
immune.part <- RunPCA(object = immune.part, npcs = 5, verbose = TRUE)
immune.part <- RunUMAP(object = immune.part, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.part <- FindNeighbors(object = immune.part, dims = 1:5, verbose = TRUE)
immune.part <- FindClusters(object = immune.part, verbose = TRUE, resolution = 0.8)
DimPlot(immune.part)
a = DoHeatmap(subset(immune.part, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '03 Stroma DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.part = subset(immune.part,seurat_clusters %in% 
                       c(0:22,25,31))

celltype.marker = c('CD31','Podoplanin',# Endothelial
                    'SMA',
                    'CD45','CD3e','CD20','CD79a','CD68','MPO' # Immune marker
)
VariableFeatures(immune.part) <- celltype.marker
immune.part = ScaleData(immune.part)
immune.part <- RunPCA(object = immune.part, npcs = 5, verbose = TRUE)
immune.part <- RunUMAP(object = immune.part, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.part <- FindNeighbors(object = immune.part, dims = 1:5, verbose = TRUE)
immune.part <- FindClusters(object = immune.part, verbose = TRUE, resolution = 0.8)
DimPlot(immune.part)
a = DoHeatmap(subset(immune.part, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '04 Stroma DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

immune.part.part1 = subset(immune.part,seurat_clusters %in% c(12,15,17,21,28))
immune.part.file.part1 = data.sub[colnames(immune.part.part1),]
immune.part.file.part1$cluster =paste0("Cluster_",immune.part.part1$seurat_clusters) 
write.csv(immune.part.file.part1,'immunepartfilepart2.csv')

immune.part.tumor = subset(immune.part,seurat_clusters %in% c(0))
immune.part.stroma = subset(immune.part,seurat_clusters %in% c(2,5,6,10,11,13,19,20,22,30,31))
immune.part.lymphoid = subset(immune.part,seurat_clusters %in% c(1,3,7,8,14,18,23,25,27,17))
immune.part.myeloid = subset(immune.part,seurat_clusters %in% c(9,16,24))


# lymphoid
immune.marker = c('CD3e','CD20','CD4','CD8','CD79a','Granzyme.B','PD.1','IgA','IgG','FOXP3')
VariableFeatures(immune.part.lymphoid) <- immune.marker
#immune.part.lymphoid@assays$RNA@scale.data[1:5,1:5]
immune.part.lymphoid = ScaleData(immune.part.lymphoid)
immune.part.lymphoid <- RunPCA(object = immune.part.lymphoid, npcs = 5, verbose = TRUE)
immune.part.lymphoid <- RunUMAP(object = immune.part.lymphoid, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.part.lymphoid <- FindNeighbors(object = immune.part.lymphoid, dims = 1:5, verbose = TRUE)
immune.part.lymphoid <- FindClusters(object = immune.part.lymphoid, verbose = TRUE, resolution = 0.5)
DimPlot(immune.part.lymphoid)
a = DoHeatmap(subset(immune.part.lymphoid, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '05 lymphoid DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

cd4t =  subset(immune.part.lymphoid,seurat_clusters %in% c(4,7,10,13,14,15))
cd8t =  subset(immune.part.lymphoid,seurat_clusters %in% c(3,6,12))
b =  subset(immune.part.lymphoid,seurat_clusters %in% c(2,9,11,17))

lymphoid.mix =  subset(immune.part.lymphoid,seurat_clusters %in% c(0,1,5,8))
immune.marker = c('CD3e','CD20','CD4','CD8','CD79a','Granzyme.B','PD.1','IgA','IgG','FOXP3')
VariableFeatures(lymphoid.mix) <- immune.marker
#immune.part.lymphoid@assays$RNA@scale.data[1:5,1:5]
lymphoid.mix = ScaleData(lymphoid.mix)
lymphoid.mix <- RunPCA(object = lymphoid.mix, npcs = 5, verbose = TRUE)
lymphoid.mix <- RunUMAP(object = lymphoid.mix, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
lymphoid.mix <- FindNeighbors(object = lymphoid.mix, dims = 1:5, verbose = TRUE)
lymphoid.mix <- FindClusters(object = lymphoid.mix, verbose = TRUE, resolution = 0.5)
DimPlot(lymphoid.mix)
a = DoHeatmap(subset(lymphoid.mix, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '06 lymphoid mix DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

lymphoid.mix.cd4t =  subset(lymphoid.mix,seurat_clusters %in% c(4,9,10))
lymphoid.mix.cd8t =  subset(lymphoid.mix,seurat_clusters %in% c(1,2,12,13))
lymphoid.mix.b =  subset(lymphoid.mix,seurat_clusters %in% c(6,7,11))

# myeloid
rownames(immune.part.myeloid)
immune.marker = c('CD14','CD163','CD68','MPO','S100A8-9')
VariableFeatures(immune.part.myeloid) <- immune.marker
immune.part.myeloid = ScaleData(immune.part.myeloid)
immune.part.myeloid <- RunPCA(object = immune.part.myeloid, npcs = 3, verbose = TRUE)
immune.part.myeloid <- RunUMAP(object = immune.part.myeloid, dims = 1:3, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.part.myeloid <- FindNeighbors(object = immune.part.myeloid, dims = 1:3, verbose = TRUE)
immune.part.myeloid <- FindClusters(object = immune.part.myeloid, verbose = TRUE, resolution = 0.5)
DimPlot(immune.part.myeloid)
a = DoHeatmap(subset(immune.part.myeloid, downsample=2000), features = immune.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '07 immune myeloid DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

macrophage = subset(immune.part.myeloid,seurat_clusters %in% c(0,1,3,8,13,15,16,17))
neutrophil = subset(immune.part.myeloid,seurat_clusters %in% c(2,4,5,6,7,9,11,12))

# stroma
rownames(immune.part.stroma)
immune.marker = celltype.marker = c('CD31','Podoplanin',# Endothelial
                                    'SMA',
                                    'CD45','CD3e','CD20','CD79a','CD68','MPO' # Immune marker
)
VariableFeatures(immune.part.stroma) <- immune.marker
#immune.part.lymphoid@assays$RNA@scale.data[1:5,1:5]
immune.part.stroma = ScaleData(immune.part.stroma)
immune.part.stroma <- RunPCA(object = immune.part.stroma, npcs = 5, verbose = TRUE)
immune.part.stroma <- RunUMAP(object = immune.part.stroma, dims = 1:5, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
immune.part.stroma <- FindNeighbors(object = immune.part.stroma, dims = 1:5, verbose = TRUE)
immune.part.stroma <- FindClusters(object = immune.part.stroma, verbose = TRUE, resolution = 0.5)
DimPlot(immune.part.stroma)

a = DoHeatmap(subset(immune.part.stroma, downsample=2000), features = c('SMA','CD31','Podoplanin'), angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '07 immune stroam DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

fibroblast = subset(immune.part.stroma,seurat_clusters %in% c(0,3,4,7,10,12,13,14,19))
endothelial = subset(immune.part.stroma,seurat_clusters %in% c(1,2,5,6,8,9,11,16,17,18,20))

# tumor
tumor.part = data.obj[,c(colnames(tumor,mix.tumor,immune.part.tumor))]
rownames(tumor.part)
celltype.marker = c('EpCAM','Pan.Cytokeratin','Hepar1')
VariableFeatures(tumor.part) <- celltype.marker
tumor.part = ScaleData(tumor.part)
tumor.part <- RunPCA(object = tumor.part, npcs = 2, verbose = TRUE)
tumor.part <- RunUMAP(object = tumor.part, dims = 1:2, verbose = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)
tumor.part <- FindNeighbors(object = tumor.part, dims = 1:2, verbose = TRUE)
tumor.part <- FindClusters(object = tumor.part, verbose = TRUE, resolution = 0.3)
DimPlot(tumor.part)
a = DoHeatmap(subset(tumor.part, downsample=2000), features = celltype.marker, angle = 90, size = 4, disp.min = -1, disp.max = 1)
pdf(file = '08 Tumor DoHeatmap marker.pdf',height = 7,width = 20)
print(a)
dev.off()

tumor.part.part = subset(tumor.part,seurat_clusters %in% c(0:39))
tumor.part.part.file = data.sub[colnames(tumor.part.part),]
tumor.part.part.file$cluster =paste0("Cluster_",tumor.part.part$seurat_clusters) 
write.csv(tumor.part.part.file,'tumorcheck3.csv')

bile.duct = subset(tumor.part,seurat_clusters %in% c(3,12,15,16,23,27,28,30,32,33,34))
tumor = tumor.part[, !colnames(tumor.part) %in% colnames(bile.duct)]

#total
bile.duct = subset(tumor.part,seurat_clusters %in% c(3,12,15,16,23,27,28,30,32,33,34))
tumor = tumor.part[, !colnames(tumor.part) %in% colnames(bile.duct)]

fibroblast = subset(immune.part.stroma,seurat_clusters %in% c(0,3,4,7,10,12,13,14,19))
endothelial = subset(immune.part.stroma,seurat_clusters %in% c(1,2,5,6,8,9,11,16,17,18,20))

macrophage = subset(immune.part.myeloid,seurat_clusters %in% c(0,1,3,8,13,15,16,17))
neutrophil = subset(immune.part.myeloid,seurat_clusters %in% c(2,4,5,6,7,9,11,12))

cd4t =  subset(immune.part.lymphoid,seurat_clusters %in% c(4,7,10,13,14,15))
cd8t =  subset(immune.part.lymphoid,seurat_clusters %in% c(3,6,12))
b =  subset(immune.part.lymphoid,seurat_clusters %in% c(2,9,11,17))

data.sub[colnames(b),'celltype'] = 'B_cell'
data.sub[colnames(cd8t),'celltype'] = 'CD8T'
data.sub[colnames(cd4t),'celltype'] = 'CD4T'
data.sub[colnames(macrophage),'celltype'] = 'Macrophage'
data.sub[colnames(neutrophil),'celltype'] = 'Neutrophil'
data.sub[colnames(endothelial),'celltype'] = 'Endothelial'
data.sub[colnames(fibroblast),'celltype'] = 'Fibroblast'
data.sub[colnames(tumor),'celltype'] = 'Tumor'
data.sub[colnames(bile.duct),'celltype'] = 'Bile_Duct'
data.sub[is.na(data.sub$celltype),'celltype'] = 'Unassign'
table(data.sub$celltype)
write.csv(data.sub,'...\\Add celltype1.csv')



filtered.cor <- readRDS('Filtered matched TLS.rds')
ta.tls.hd <- readRDS('CODEX-identified ta-TLS group.rds')
ta.tls.hd <- subset(ta.tls.hd,group == 'T-TLS')
common.tls <- intersect(rownames(filtered.cor),rownames(ta.tls.hd))

count.total <- readRDS('HD Count.rdata')
expr <- CreateSeuratObject(count.total)
expr <- NormalizeData(expr)
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr,  verbose = F)
PCS = 10
k.param = 10
expr <- RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr <- FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr <- RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr <- FindClusters(expr,resolution = 0.3,verbose = FALSE)
DimPlot(expr,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
expr <- expr[,common.tls]
DotPlot(expr,features = c('CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','CXCL7','CXCL8',
                          'CXCL12','CCL3','CCL4','CCL5'))
expr$orig.ident
expr$group = 'one'
avg_expr <- AverageExpression(expr, features = c('CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','CXCL8',
                                                 'CXCL12','CCL3','CCL4','CCL5','CCL7','CCL8','CCL15','CCL23'),
                              group.by = 'group')
avg_expr <- AverageExpression(expr, features = c('CCL3','CCL5','CCL7','CCL8','CCL15','CCL23'),
                              group.by = 'group')
a <- pheatmap(as.matrix(avg_expr$RNA),cellwidth = 20,cellheight = 10,
             cluster_rows = FALSE, cluster_cols = FALSE,
             color = colorRampPalette(c("lightgrey", "#7971ea"))(100),
             main = "Chemoattractant Expression in Each TLS")


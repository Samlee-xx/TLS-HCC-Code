
###### Pseudotime analysis of HD TLS ######
library(monocle)

setwd('...\\06 HD TLS Pseudotime Analysis')
colnames(hd.135121.expr) = paste0('135121_',colnames(hd.135121.expr))
colnames(hd.164697.expr) = paste0('164697_',colnames(hd.164697.expr))
colnames(hd.178031.expr) = paste0('178031_',colnames(hd.178031.expr))
colnames(hd.238966.expr) = paste0('238966_',colnames(hd.238966.expr))
colnames(hd.239117.expr) = paste0('239117_',colnames(hd.239117.expr))

common.genes <- Reduce(intersect, list(
  rownames(hd.135121.expr),
  rownames(hd.164697.expr),
  rownames(hd.178031.expr),
  rownames(hd.238966.expr),
  rownames(hd.239117.expr)
))

hd.135121.expr = hd.135121.expr[common.genes,]
hd.164697.expr = hd.164697.expr[common.genes,]
hd.178031.expr = hd.178031.expr[common.genes,]
hd.238966.expr = hd.238966.expr[common.genes,]
hd.239117.expr = hd.239117.expr[common.genes,]

count.total = cbind(hd.135121.expr,hd.164697.expr,hd.178031.expr,hd.238966.expr,hd.239117.expr)
saveRDS(count.total,
        file = '...\\HD Count.rdata')

count.total = readRDS('...\\HD Count.rdata')
col_names = colnames(count.total)
id_numbers <- str_extract(col_names, "^\\d+")
id_letters <- str_extract(col_names, "(?<=_)[A-Z]+")

cell_metadata=data.frame(cell=colnames(count.total),
                         location = id_letters,
                         Sample = id_numbers)
cell_metadata$stage = case_when(cell_metadata$Sample %in% c('164697','178031')~'Early',
                                cell_metadata$Sample %in% c('135121','176125')~'Middle',
                                cell_metadata$Sample %in% c('238966','239117')~'Late')

rownames(cell_metadata)=cell_metadata$cell
gene_metadata=data.frame(gene=rownames(count.total))
rownames(gene_metadata)=gene_metadata$gene
gene_metadata$gene_short_name=gene_metadata$gene

pd <- new('AnnotatedDataFrame', data = cell_metadata) 
fd <- new('AnnotatedDataFrame', data = gene_metadata)

cds <- newCellDataSet(as.matrix(count.total),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "stage")  + 
  theme(legend.key = element_rect(fill=NA),
        legend.title = element_text(color = 'chocolate',
                                    size=14,face=2))+
  guides(color=guide_legend(override.aes = list(size=6)))

a = plot_cell_trajectory(cds,color_by = "Pseudotime")
saveRDS(cds,'HD data Pseudotime.rds')


diff_test_res <- differentialGeneTest(cds, 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)")
sig_genes <- subset(diff_test_res, qval < 0.05)
sig_gene_names <- rownames(sig_genes)
pseudotime_values <- pData(cds)$Pseudotime
spearman_results <- apply(exprs(cds)[sig_gene_names, ], 1, function(gene_expr) {
  cor(gene_expr, pseudotime_values, method = "spearman", use = "complete.obs")
})

correlation_df <- data.frame(
  Gene = names(spearman_results),
  Spearman_Correlation = spearman_results
)
correlation_df <- correlation_df[order(correlation_df$Spearman_Correlation, decreasing = TRUE), ]
head(correlation_df, 10)

a = plot_genes_in_pseudotime(cds["", ], color_by = "stage")

gene_to_cluster <- rownames(sig_genes[order(sig_genes$qval), ])[1:100]  
gene_to_cluster = c('CD274','PDCD1','CTLA4','IDO1','FOXP3',
                    'LGALS3','TIGIT','LGALS1','LAG3','S100A8','CD177','CSF3R'
                 
)
a = plot_pseudotime_heatmap(cds[gene_to_cluster, ],
                            show_rownames = TRUE,
                            cluster_rows = FALSE,
                            cores = 4, 
                            return_heatmap = TRUE,
                            hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))

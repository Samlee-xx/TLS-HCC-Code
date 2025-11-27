
#### Cholangiocyte along the tumor from HD data ####

library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

localdir  <-  '...\\outs'
hd.data <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

bile.duct  <-  read.csv('\\sample 6 Bile Duct Barcode.csv')
colnames(hd.data)

hd.bile  <-  hd.data[,bile.duct$Barcode]
dim(hd.bile)
dim(hd.data)

barcode.xy <- read.csv('sample 6 Barcode centroids.csv')
rownames(barcode.xy) <- barcode.xy$barcode
barcode.xy <- barcode.xy[bile.duct$Barcode,]

ggplot(barcode.xy)+
  geom_point(aes(x,y),size=0.5)+
  theme_classic()

expr_matrix <- hd.bile@assays$Spatial.008um@layers$counts %>% as.matrix()
colnames(expr_matrix) = colnames(hd.bile)
rownames(expr_matrix) = rownames(hd.bile)

dim(expr_matrix)
rownames(expr_matrix)
sum(!colnames(hd.bile) == barcode.xy$barcode)
barcode_data <- data.frame(Barcode = colnames(expr_matrix), y = barcode.xy$y)

save(hd.bile, barcode.xy, file = "sample_6_Bile_Duct_HD_data.RData")
expr_data <- t(expr_matrix)
expr_data <- cbind(barcode_data, expr_data)

ggscatter(expr_data, x = 'x', y = 'MMP7',
          fill ='#8491B4FF',color='darkgrey',size = 1,shape = 21, # Points color, shape and size
          add = "reg.line",# Add regressin line, # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = T, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x.npc = "left", 
                                label.y.npc = "top",label.sep = "\n"))

ggplot(expr_data)+
  geom_point(aes(x,MMP7),size=0.5)+
  theme_minimal()

pdf(file = 'E:\\01-TLS\\03-HD\\05-Analysis\\03-AlongAxis\\MMP7.pdf',
    height = 3,width = 3)
print(a)
dev.off()

# calculate each gene along x or y axis
gene_correlations <- apply(expr_matrix, 1, function(gene_expr) {
  cor.test(barcode.xy$x, gene_expr, method = "spearman")
})

cor_results <- data.frame(
  Gene = rownames(expr_matrix),
  Correlation = sapply(gene_correlations, function(x) x$estimate),
  PValue = sapply(gene_correlations, function(x) x$p.value)
)

positive_percentage <- apply(expr_matrix, 1, function(gene_expr) {
  sum(gene_expr > 0) / length(gene_expr) * 100
})

cor_results$PositivePercentage <- positive_percentage

significant_genes <- cor_results %>%
  filter(Correlation > 0, PValue < 0.05) %>%
  arrange(desc(Correlation))


library(clusterProfiler)
library(org.Hs.eg.db)

significant_genes <- subset(significant_genes,PValue <= 0.05)

gene.id <- mapIds(org.Hs.eg.db, 
                  keys = significant_genes$Gene,keytype = "SYMBOL", column="ENTREZID")
GO <- enrichGO(gene = gene.id,
            OrgDb = org.Hs.eg.db,
            pvalueCutoff =0.05,	
            qvalueCutoff = 0.05,	
            ont="BP",	
            readable =T)

Go.sub <- GO@result
Go.sub <- filter(Go.sub,p.adjust <=0.05)
Go.sub <- subset(Go.sub,Description %in% c('response to chemokine',
                                          'neutrophil migration',
                                          'regulation of neutrophil chemotaxisy',
                                          'neutrophil chemotaxis',
                                          'granulocyte chemotaxis',
                                          'granulocyte migration'
))

b <- lapply(str_split(Go.sub$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )

Go.sub$GeneRatio = b

mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
a <- ggplot(data =Go.sub, aes(x =GeneRatio , y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio", y = "", title = "") +
  geom_text(aes(x = 0.000, 
                label = Description),
            hjust = 0)+
  theme_classic() + mytheme2

pdf(file = 'E:\\01-TLS\\03-HD\\05-Analysis\\03-AlongAxis\\238966 Bile Duct along axis go terms2.pdf',height = 3,width = 6)
print(a)
dev.off()

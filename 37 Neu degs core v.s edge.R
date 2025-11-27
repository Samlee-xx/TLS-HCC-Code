

neu.hcc <- readRDS('01 Neu seurat.rds')
unique(neu.hcc$celltype)
marker <- FindMarkers(neu.hcc,ident.1 = c('Neu_S2_MIF'),ident.2 = c('Neu_S0_S100A8','Neu_S4_S100A12','Neu_S1_APOC3'))
marker$gene <- rownames(marker)

df <- marker
df$logpvalue <- -log10(df$p_val_adj)
df$group <- "none"
df$group[which((df$avg_log2FC >=0.5) & (df$p_val_adj < 0.05))] <- "Up"  
df$group[which((df$avg_log2FC <= -0.5) & (df$p_val_adj < 0.05))] <- "Down" 
df <- df[complete.cases(df$gene),]
df$label  <-  ''
gene <- c('CXCR2','CCR1','MIF','LGALS3','CD74','CCL3','CCL4','CCRL2',
         'S100A8','S100A9','S100A12','HLA-DRA')
df$label[match(gene,df$gene)] <- gene 
df$color <- ifelse(df$group == "none" & df$label == "", "color1",   
                   ifelse(df$group == "Up" & df$label == "", "color2",   
                          ifelse(df$group == "Down" & df$label == "", "color3", 
                                 ifelse(df$group == "Up" & df$label != "", "color4", "color5"))))  
df <- arrange(df, color)  
colnames(df)
library(ggpubr)
a <- ggscatter(df,
              x="avg_log2FC", 
              y="logpvalue", 
              color = "color",  
              palette = c("#bcbcbc","#ffab84","#8abddc","#be000e","#0051a6"), 
              label = df$label,   
              font.label = c(15,"plain","black"),  
              repel = T ) +  
  labs(title="Neutrophil Chemokines")+   
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+  
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        legend.position = 'none')  


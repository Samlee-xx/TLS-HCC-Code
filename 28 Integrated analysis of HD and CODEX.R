
###### Celltypes and function correlation analysis ######

library(GO.db)
library(org.Hs.eg.db)
library(msigdbr)

all_gene_sets <- msigdbr_collections() 
print(unique(all_gene_sets$gs_cat))

filtered.cor = readRDS('...\\Filtered matched TLS.rds')
ta.tls.hd = readRDS('...\\CODEX-identified ta-TLS group.rds')
ta.tls.hd = subset(ta.tls.hd,group == 'T-TLS')
common.tls = intersect(rownames(filtered.cor),rownames(ta.tls.hd))
ta.tls.hd.meta = ta.tls.hd[common.tls,]

combat_corrected_data = readRDS('...\\TLS HD.rds')
colnames(combat_corrected_data)
combat_corrected_data = combat_corrected_data[,rownames(ta.tls.hd.meta)]

sign <- read.delim("...\\04 Immune Signature.csv",sep=',')
colnames(sign)
gene.list = list()
for (i in c(1,2,3,4,5,6,8,12,15,16,17,18,19,20,51,52,53,54,55,56,57,58,59,60,61,62,72,
            73,74,75,76,77,78,79,90,91,92,93,94,95,96,97)){
  gene.list[[colnames(sign)[i]]] = unique(sign[,i])
}

genesets <- msigdbr(
  species = "Homo sapiens",     
  category = "C5",                   
  #subcategory = "CP:GO"            
) %>% 
  filter(grepl("GOBP_NEUTROPHIL_CHEMOTAXIS", gs_name, ignore.case = TRUE)) 

gene.list[['Neutrophil_Chemokines']] = c('CXCL1','CXCL2','CXCL3','CXCL6','CXCL8','CXCL5')
sig_tme <- calculate_sig_score(pdata           = NULL,
                               eset            = combat_corrected_data,
                               signature       = gene.list,
                               method          = "ssgsea",
                               mini_gene_count = 0,
                               adjust_eset = TRUE)

sig_tme = as.data.frame(sig_tme)
rownames(sig_tme) = sig_tme$ID
sig_tme = sig_tme[,-c(1,2)]
rownames(sig_tme) == rownames(ta.tls.hd.meta)

setwd('...\\0306 second')
plot.table = cbind(sig_tme,ta.tls.hd.meta)

colnames(plot.table)
heat.plot = plot.table[,c(1,3,4,5,6,10,11,14,20,22,23,24,25,26,28,39,40,42,43:51)]
colnames(heat.plot)[9] = 'T_receptor_signaling'
colnames(heat.plot)[10] = 'CD40 pathway'
colnames(heat.plot)[11] = 'CTL pathway'
colnames(heat.plot)[12] = 'IL2 pathway'
colnames(heat.plot)[13] = 'IL4 pathway'
colnames(heat.plot)[14] = 'IL10 pathway'
colnames(heat.plot)[15] = 'GC B cell'
colnames(heat.plot)[16] = 'B activation'
colnames(heat.plot)[17] = 'B cell receptor signaling'
colnames(heat.plot)[18] = 'Lymphocyte costimulation'
colnames(heat.plot)

heat.plot = heat.plot[,c(2,3,5,6,7,9,10,11,12,13,14,15,16,17,18:27)]
colnames(heat.plot)

pathway_data <- heat.plot[, c(1,7,8,13,14,16)]
celltype_data <- heat.plot[, 17:ncol(heat.plot)] 
cor_matrix <- cor(celltype_data, pathway_data, method = "spearman")

a=pheatmap(
  cor_matrix,
  color = colorRampPalette(c("#739559", "white", '#f50538'))(50),  
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  fontsize_row = 10, 
  fontsize_col = 10, 
  main = "Spearman Correlation: Pathways vs. Cell Proportions"
)






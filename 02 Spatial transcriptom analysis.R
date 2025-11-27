
###### This code used for spatial transcriptom analysis ----

## 1.TLS Subtypes and cell2location  ----

# TLS Subtypes

tls=readRDS('...\\tls.rds') 

a = DimPlot(tls)+theme(axis.line=element_line(colour='black',size=0.5,lineend = 'square'))+
  guides(color=guide_legend(override.aes = list(size=6)))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        legend.direction = "vertical",
        legend.position = "right",
        #legend.box = "vertical"
  )+ggtitle('TLS Subtypes')+
  scale_color_manual(values =c('#EA966B',"#DA8DC0", '#AFD766'))

pdf(file = '...\\01 tls_umap.pdf',height = 3,width = 4.5)
print(a)
dev.off()

# cell2location run by python, now we analysis the result 


setwd('...\\06 Cell2location')
data = read.csv('tls_cell2location_out.csv')
rownames(data) = data$X
data = data[,c(6:16)]
data = data[,-c(2,7)]

data$b.linkage = data$B.Cell+data$Plasma.Cell
data$t.linkage = data$CD4.T.Cell+data$CD8.T.Cell
colnames(data)
data = data[,-c(1,2,3,9)]
colnames(data)

normalize_to1 = function(x){
  x/sum(x)*100
}
data.norm = apply(data,1,normalize_to1)
data.norm =  data.norm %>% t() %>% as.data.frame()

sum(!colnames(tls) == rownames(data.norm)) #0
data.norm$group = tls$seurat_clusters
data.norm = melt(data.norm)

unique(data.norm$variable)
celltypes = unique(data.norm$variable)
data.norm.cell = subset(data.norm,variable == celltypes[7])

a = ggboxplot(data.norm.cell, x = "group", y = "value", fill = "group",
              ylab = "Value",
              title = celltypes[1],
              outlier.shape = NA)+
  stat_compare_means(method = "anova", label = "p.format")+
  #geom_jitter(width = 0.3, size = 1)+
  geom_hline(yintercept = median(data.norm.cell$value),color='red',linetype='dashed')+
  #ylim(0,25)+
  scale_fill_manual(values =  c('#FF9900','#FF0033','#339933'))

pdf(file = 'T annova.pdf',height = 2.7,width =3.5)
print(a)
dev.off()

## 2.TLS DEGs using DESeq2 ----
### T-TLS
tls.count = readRDS('...\\TLS count matrix.rds')
tls.group <- ifelse(tls$seurat_clusters == "T-enriched", "T-TLS", "others")
tls.group <- factor(tls.group, levels = c("others", "T-TLS"))
sum(!colnames(tls.count) == colnames(tls)) #0
colData <- data.frame(group = tls.group)
rownames(colData) <- colnames(tls.count)

dds <- DESeqDataSetFromMatrix(countData = tls.count,
                              colData = colData,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "T-TLS", "others"))
res <- res[order(res$pvalue), ] %>% 
  na.omit() %>%
  as.data.frame()
res$gene <- rownames(res)
t.tls.up <- res %>%
  filter(padj <= 0.01, log2FoldChange > 0)

write.table(t.tls.up,
            file = "...\\T_TLS_Gene_Sig.csv",
            sep = ",",
            row.names = FALSE,  
            quote = FALSE)   

t.tls.up$enterz <- mapIds(org.Hs.eg.db, keys = t.tls.up$gene, keytype = "SYMBOL", column="ENTREZID")
GO=enrichGO(gene = t.tls.up$enterz,
            OrgDb = org.Hs.eg.db, 
            pvalueCutoff =0.05,
            qvalueCutoff = 0.05,	
            ont="BP",	
            readable =T)	
t.go.result = GO@result
t.go.result = filter(t.go.result,p.adjust <=0.05)

t.go.result$GeneRatio <- paste0("'", t.go.result$GeneRatio) 
write.table(t.go.result, "...\\T_TLS_go_terms.csv", sep = ",", row.names = FALSE, quote = TRUE)


go.filtered <- t.go.result %>%
  filter(grepl("T cell|immune|chemotaxis|antigen|neutrophil|cytotoxicity|cytokine", Description, ignore.case = TRUE))
a = enrichmentNetwork(go.filtered, 
                      colorBy = 'pvalue', 
                      colorType = 'pval', nodeSize = "Count",drawEllipses = TRUE)+
  scale_color_gradientn(colours = c("#B83D3D",'white','#1A5592'),
                        name = "pvalue")

pdf(file='...\\T-TLS GO plot.pdf',height = 7,width =7)
print(a)
dev.off()

### B-TLS
tls.group <- ifelse(tls$seurat_clusters == "B-enriched", "B-TLS", "others")
tls.group <- factor(tls.group, levels = c("others", "B-TLS"))
all(colnames(tls.count) == colnames(tls))

colData <- data.frame(group = tls.group)
rownames(colData) <- colnames(tls.count)
dds <- DESeqDataSetFromMatrix(countData = tls.count,
                              colData = colData,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "B-TLS", "others")) %>%
  as.data.frame() %>%
  na.omit() %>%
  arrange(pvalue)
res$gene <- rownames(res)
b.tls.up <- res %>%
  filter(padj <= 0.01, log2FoldChange > 0)

write.table(b.tls.up,
            file = "...\\B_TLS_Gene_Sig.csv",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

b.tls.up$enterz <- mapIds(org.Hs.eg.db,
                          keys = b.tls.up$gene,
                          keytype = "SYMBOL",
                          column = "ENTREZID",
                          multiVals = "first")

GO <- enrichGO(gene = na.omit(b.tls.up$enterz),
               OrgDb = org.Hs.eg.db,
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               ont = "BP",
               readable = TRUE)

b.go.result <- GO@result %>%
  filter(p.adjust <= 0.05)

b.go.result$GeneRatio <- paste0("'", b.go.result$GeneRatio) 
write.table(b.go.result, "...\\B_TLS_go_terms.csv", sep = ",", row.names = FALSE, quote = TRUE)


a = enrichmentNetwork(b.go.result, 
                      colorBy = 'pvalue', 
                      colorType = 'pval', nodeSize = "Count",drawEllipses = TRUE)+
  scale_color_gradientn(colours = c("#B83D3D",'white','#1A5592'),
                        name = "pvalue")

pdf(file='...\\B-TLS GO plot.pdf',height =4,width =4)
print(a)
dev.off()

## Metabolism subtype
tls.count = readRDS('...\\TLS count matrix.rds')
tls.group <- ifelse(tls$seurat_clusters == "M-TLS", "M-TLS", "others")
table(tls$seurat_clusters)
tls.group <- factor(tls.group, levels = c("others", "M-TLS"))
sum(!colnames(tls.count) == colnames(tls)) #0
colData <- data.frame(group = tls.group)
rownames(colData) <- colnames(tls.count)

dds <- DESeqDataSetFromMatrix(countData = tls.count,
                              colData = colData,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "M-TLS", "others"))
res <- res[order(res$pvalue), ] %>% 
  na.omit() %>%
  as.data.frame()
res$gene <- rownames(res)

m.tls.up <- res %>%
  filter(padj <= 0.01, log2FoldChange > 0)

write.table(m.tls.up,
            file = "...\\M_TLS_Gene_Sig.csv",
            sep = ",",
            row.names = FALSE,   
            quote = FALSE)   

m.tls.up$enterz <- mapIds(org.Hs.eg.db, keys = m.tls.up$gene, keytype = "SYMBOL", column="ENTREZID")
GO=enrichGO(gene = m.tls.up$enterz,
            OrgDb = org.Hs.eg.db, 
            pvalueCutoff =0.05,
            qvalueCutoff = 0.05,	
            ont="BP",	
            readable =T)	
m.go.result = GO@result
m.go.result = filter(m.go.result,p.adjust <=0.05)

m.go.result$GeneRatio <- paste0("'", m.go.result$GeneRatio) 
write.table(m.go.result, "...\\M_TLS_go_terms.csv", sep = ",", row.names = FALSE, quote = TRUE)

a = enrichmentNetwork(m.go.result, 
                      colorBy = 'pvalue', 
                      colorType = 'pval', nodeSize = "Count",drawEllipses = TRUE)+
  scale_color_gradientn(colours = c("#B83D3D",'white','#1A5592'),
                        name = "pvalue")

pdf(file='...\\M-TLS GO plot.pdf',height = 7,width =7)
print(a)
dev.off()

## Save gene signature 
saveRDS(b.tls.up,'E:\\01-TLS\\01-ST\\B_TLS_Sig.rds')
saveRDS(t.tls.up,'E:\\01-TLS\\01-ST\\T_TLS_Sig.rds')
saveRDS(m.tls.up,'E:\\01-TLS\\01-ST\\M_TLS_Sig.rds')

b.tls.up = readRDS('E:\\01-TLS\\01-ST\\B_TLS_Sig.rds')
t.tls.up = readRDS('E:\\01-TLS\\01-ST\\T_TLS_Sig.rds')

## 3.TLS Prognosis and immunotherapy ----
library(IOBR)
library(dplyr)
library(reshape2)
library(survival)
library(survminer)
### 3.1 Bulk data from TCGA
## calculate tls score
total.count <- readRDS("...\\total_tpm_408.rds")
rownames(total.count) <- total.count$gene_name
total.count <- total.count[, -1]

gene.total <- list(TLS.0 = t.tls.up$gene, TLS.1 = b.tls.up$gene)

sig_tme <- calculate_sig_score(
  pdata           = NULL,
  eset            = total.count,
  signature       = gene.total,
  method          = "ssgsea",
  mini_gene_count = 0,
  adjust_eset     = TRUE
)

colnames(sig_tme)[1] <- "Sample.ID"
## make up clinical data
sample.sheet <- read.table('...\\gdc_sample_sheet.2022-12-03.tsv',
                           sep = '\t', quote = "", header = TRUE)
info.patient <- read.table('...\\clinical.cart.2022-12-02\\clinical.tsv',
                           sep = '\t', quote = "", header = TRUE)

sample.sheet <- sample.sheet %>%
  dplyr::select(File.Name, Sample.ID, Sample.Type, Case.ID)

info.patient <- info.patient %>%
  dplyr::select(case_submitter_id, days_to_death, vital_status, ajcc_pathologic_stage) %>%
  rename(Case.ID = case_submitter_id)

info <- right_join(sample.sheet, info.patient, by = "Case.ID")
final <- left_join(sig_tme, info, by = "Sample.ID") %>%
  distinct(Sample.ID, .keep_all = TRUE)

## clean up clinical data
final$days_to_death <- gsub("'--", NA, final$days_to_death)
final$ajcc_pathologic_stage <- gsub("'--", NA, final$ajcc_pathologic_stage)

final$vital_status <- gsub("Not Reported", NA, final$vital_status)
final$vital_status <- gsub("Alive", 0, final$vital_status)
final$vital_status <- gsub("Dead", 1, final$vital_status)

final$days_to_death <- as.numeric(final$days_to_death)
final$vital_status <- as.numeric(final$vital_status)

final <- final %>%
  filter(!is.na(vital_status), Sample.Type == "Primary Tumor")

## prognosis analysis
table(final$ajcc_pathologic_stage)
final.early <- final %>%
  filter(ajcc_pathologic_stage %in% c("Stage I", "Stage II"))

sur.cut.early <- surv_cutpoint(final.early,
                               time = "days_to_death",
                               event = "vital_status",
                               variables = "TLS.0",
                               minprop = 0.45)
summary(sur.cut.early)

sur.cat.early <- surv_categorize(sur.cut.early)
names(sur.cat.early) <- c("OSday", "Death", "group")

fit.early <- survfit(Surv(OSday, Death) ~ group, data = sur.cat.early)

a = ggsurvplot(fit.early,
               pval = TRUE, conf.int = FALSE,
               risk.table.col = "strata",
               risk.table = TRUE,
               linetype = "strata",
               surv.median.line = "hv",
               ggtheme = theme_classic(),
               palette = c('#d62728', '#1f77b4')) +
  ggtitle("T-enriched TLS (Early Stage: I/II)")

sur.cat.early$group <- factor(sur.cat.early$group, levels = c("low", "high"))
cox.early <- coxph(Surv(OSday, Death) ~ group, data = sur.cat.early)
summary(cox.early)

pdf(file='...\\B-enriched TLS.pdf',height = 5,width =4)
print(a)
dev.off()

### 3.2 ICB from DongChen Cancer cell

setwd('...\\GSE235863')
bulk.data <- read.table('GSE235863_bulk_rna_seq_tpm.txt',
                        sep = '\t', header = TRUE, row.names = 1)

bulk.data <- bulk.data[!duplicated(bulk.data$geneName), ]
rownames(bulk.data) <- bulk.data$geneName
bulk.data <- bulk.data[, -1]

gene.total <- list(TLS.0 = t.tls.up$gene,
                   TLS.1 = b.tls.up$gene)
sig_tme <- calculate_sig_score(pdata           = NULL,
                               eset            = bulk.data,
                               signature       = gene.total,
                               method          = "zscore",
                               mini_gene_count = 0,
                               adjust_eset     = TRUE)

sig_tme$group <- c('NR','PR','PR','NR','CR','NR','PR','PR','PR','NR','CR','PR','CR','CR','CR')
comparisons <- list(c("NR", "PR"), c("NR", "CR"), c("PR", "CR"))


x <- "group"
y <- "TLS.1"

stat.test <- compare_means(
  formula = as.formula(paste(y, "~", x)),
  data = sig_tme,
  method = "wilcox.test",
  p.adjust.method = "holm",
  comparisons = comparisons
)

# set position for sig format
ymax <- max(sig_tme[[y]], na.rm = TRUE)
stat.test$y.position <- seq(ymax * 1.05,
                            by = ymax * 0.05,
                            length.out = nrow(stat.test))

colnames(sig_tme)
p_tls0 <- ggboxplot(sig_tme, x = 'group', y = 'TLS.1',
                    color = "group", palette = "jco", add = "jitter") +
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01) +
  theme_classic() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("B-enriched TLS")

pdf("B_enriched_TLS_adjusted.pdf", height = 4, width = 3)
print(p_tls0)
dev.off()

### 3.3 Imbrave150
#### CR PR SD PD
setwd('...\\03 Imbrave50')

pheno.data <- read.table('...\\IMbrave150_GO30140\\IMbrave150.cli.txt',
                         sep='\t',
                         header = TRUE)

table(pheno.data$Treatment)
table(pheno.data$Treatment,pheno.data$Rsponse)

expr.data <-  read.table('...\\IMbrave150_GO30140\\IMbrave150.exp.txt',
                         sep ='\t',header = TRUE)
rownames(expr.data) = expr.data$Symbol
expr.data = expr.data[,-c(1,2)] 

pheno.data$anon_sampleId = gsub('-','.',pheno.data$anon_sampleId)
expr.data = expr.data[,pheno.data$anon_sampleId]
sum(!colnames(expr.data) %in% pheno.data$anon_sampleId)

# Signature
b.tls.up = readRDS('...\\B_TLS_Sig.rds')
t.tls.up = readRDS('...\\T_TLS_Sig.rds')

gene.total <- list(TLS.0 = t.tls.up$gene, TLS.1 = b.tls.up$gene)
sig_tme <- calculate_sig_score(pdata = NULL,
                               eset = expr.data,
                               signature = gene.total,
                               method = "zscore",
                               mini_gene_count = 0,
                               adjust_eset = TRUE)

colnames(pheno.data)[3] <- 'sample.id'
colnames(sig_tme)[2] <- 'sample.id'
sig_tme <- sig_tme %>% left_join(pheno.data, by = "sample.id")
sig_tme[sig_tme == ""] <- NA

# filtered patients with clinical outcome
sig.tme.ab <- sig_tme %>%
  filter(Treatment == 'Atezolizumab+Bevacizumab',
         Confirmed.Response_IRF %in% c('CR', 'PR', 'SD', 'PD'))

my_comparisons <- list(c("CR", "PD"), c("CR", "PR"), c("CR", "SD"),
                       c("PD", "PR"), c("PD", "SD"), c("PR", "SD"))

plot_box_with_adjusted_p <- function(df, x, y, comparisons, title = NULL, ylab = "Value") {
  stat.test <- compare_means(
    formula = as.formula(paste(y, "~", x)),
    data = df,
    method = "wilcox.test",
    p.adjust.method = "holm",
    comparisons = comparisons
  )
  stat.test$p.signif <- cut(stat.test$p.adj,
                            breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                            labels = c("****", "***", "**", "*", "ns"))
  ymax <- max(df[[y]], na.rm = TRUE)
  stat.test$y.position <- seq(ymax * 1.05, by = ymax * 0.05, length.out = nrow(stat.test))
  p <- ggboxplot(df, x = x, y = y, fill = x,
                 outlier.shape = NA, ylab = ylab, title = title) +
    geom_jitter(width = 0.3, size = 1) +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")+
    scale_fill_manual(values = c("#DA8DC0", '#AFD766', "#92A0CD",'#EA966B'))
  return(p)
}

sig.tme.ab$Confirmed.Response_IRF = factor(sig.tme.ab$Confirmed.Response_IRF,
                                           levels = c('PD','SD','PR','CR'))
a= plot_box_with_adjusted_p(df = sig.tme.ab, x = "Confirmed.Response_IRF", y = "TLS.1",
                            comparisons = my_comparisons, title = "B TLS", ylab = "Value")

a = plot_box_with_adjusted_p(df = sig.tme.ab, x = "Confirmed.Response_IRF", y = "TLS.0",
                             comparisons = my_comparisons, title = "T TLS", ylab = "Value")

pdf("T_enriched_TLS_adjusted.pdf", height = 5, width = 3)
print(a)
dev.off()


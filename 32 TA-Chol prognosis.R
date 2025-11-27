
###### TA-Chol prognosis ######
sample.6 = subset(slide.6,PValue <= 0.01)
sample.5 = subset(slide.5,PValue <= 0.01)
ta.bd.gene = union(sample.1$Gene,sample.2$Gene)
ta.bd.gene = intersect(slide.6$Gene,slide.5$Gene)

library(IOBR)
hccdb=read.table('...\\HCCDB18_mRNA_level3.txt',
                 header = TRUE,sep = '\t')
rownames(hccdb)=hccdb$Symbol
hccdb=hccdb[,-c(1,2)]
t.tls.up = readRDS('...\\T_TLS_Sig.rds')
gene.total=list(ta.bd.gene = ta.bd.gene,t.tls =t.tls.up$gene )
sig_tme <- calculate_sig_score(pdata           = NULL,
                               eset            = hccdb,
                               signature       = gene.total,
                               method          = "ssgsea",
                               mini_gene_count = 0,
                               adjust_eset = TRUE)

info1 = read.table('...\\HCCDB18.patient.txt',
                   sep='\t',header = TRUE)
rownames(info1) =info1$PATIENT_ID
info1 = as.data.frame(t(info1))
info1 = info1[-1,]
info2 = read.table('...\\HCCDB18.sample.txt',
                   sep='\t',header = TRUE)
rownames(info2) = info2$SAMPLE_ID
info2= as.data.frame(t(info2))
info2 = info2[-1,]
info2$PATIENT_ID=substr(info2$PATIENT_ID,10,15)
info1$PATIENT_ID=rownames(info1)
info1 = info1[,c("SUR","STATUS","OV_STATUS1","PATIENT_ID","TNM_STAGE_T1")]
info1$PATIENT_ID=substr(info1$PATIENT_ID,10,15)
info=info2%>% left_join(info1)
info$SAMPLE_NAME=rownames(info2)
colnames(sig_tme)[1] = 'SAMPLE_NAME' 

final=sig_tme %>% left_join(info) 
final$STATUS=gsub("Alive",0,final$STATUS)
final$STATUS=gsub("Dead",1, final$STATUS)
final$STATUS=as.numeric(final$STATUS)
final$SUR=as.numeric(final$SUR)

final$stage=case_when(final$TNM_STAGE_T1 %in% c('I','II')~'Early',
                      final$TNM_STAGE_T1 %in% c('III','IV')~'Advanced')

library(survival)
library(survminer)

table(final$TYPE)
final2 = filter(final,TYPE=='HCC')
final2 = filter(final2,stage == 'Early')

TLS50_median <- median(final2$t.tls, na.rm = TRUE)
final2$TLS50_group <- ifelse(final2$t.tls > TLS50_median, "TLS50_high", "TLS50_low")

tls_high_data <- final2 %>% filter(TLS50_group == "TLS50_high")
ta_median <- median(tls_high_data$ta.bd.gene, na.rm = TRUE)
tls_high_data$ta_group <- ifelse(tls_high_data$ta.bd.gene > ta_median, "ta_high", "ta_low")

final2$group3 <- "TLS50_low"
final2$group3[final2$TLS50_group == "TLS50_high" & final2$PATIENT_ID %in% tls_high_data$PATIENT_ID & tls_high_data$ta_group == "ta_low"] <- "TLS50_high__ta_low"
final2$group3[final2$TLS50_group == "TLS50_high" & final2$PATIENT_ID %in% tls_high_data$PATIENT_ID & tls_high_data$ta_group == "ta_high"] <- "TLS50_high__ta_high"

fit <- survfit(Surv(SUR, STATUS) ~ group3, data = final2)

a = ggsurvplot(fit,
               palette = c("#2c6e49", "#d9042b", "#262d79"),
               risk.table = TRUE,
               pval = TRUE,
               conf.int = FALSE,
               xlab = "Time in Days",
               ggtheme = theme_classic(),
               title = "Three-group Survival Analysis")

pairwise_result <- pairwise_survdiff(
  Surv(SUR, STATUS) ~ group3,
  data = final2,
  p.adjust.method = "none"  
)
pairwise_result$p.value 


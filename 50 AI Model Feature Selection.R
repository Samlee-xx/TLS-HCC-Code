
########## Feature Selection #############

library(ranger)
library(pROC)
library(randomForest)
library(ggplot2)
library(caret)
library(reshape2)
library(dplyr)

set.seed(777)
cell.meta <- readRDS('...\\TA_TLS_Metadata.rds')
cell.meta <- dcast(cell.meta,tls_id~celltype) %>% as.data.frame()
rownames(cell.meta) <- cell.meta$tls_id
cell.meta <- cell.meta[,-1]
norm=function(x){
  ((x/sum(x))*100)
} 

ct.tls <- apply(cell.meta,1,norm)
ct.tls <- t(ct.tls) %>% as.data.frame()
ct.tls$tls_id <- rownames(ct.tls)

ta.tls <- readRDS('...\\TA TLS Clustering.rds')
colnames(ta.tls)[5] = 'tls_id'
ta.tls <- ta.tls[,c('group','tls_id')]

ct.tls <- ct.tls %>% left_join(ta.tls)
ct.tls <- ct.tls %>% as.data.frame()
rownames(ct.tls) <- ct.tls$tls_id
ct.tls <- ct.tls[,-14]
colnames(ct.tls)

rm(cell.meta,ta.tls)
ta.tls.hd = readRDS('...\\CODEX-identified ta-TLS group.rds')
#merge
ct.tls$tcell = ct.tls$CD4T+ct.tls$CD8T
ta.tls.hd$tcell = ta.tls.hd$CD4T+ta.tls.hd$CD8T

colnames(ct.tls)
colnames(ta.tls.hd)

used.feature = c('B','Bile_Duct','Endothelial','Macrophage','Neutrophil','Fibroblast',
                 'tcell','group')

cohort.1 = ct.tls[,used.feature]
cohort.2 = ta.tls.hd[,used.feature]

write.csv(cohort.1,'...\\Cohort1.csv')
write.csv(cohort.2,'...\\Cohort2.csv')

cohort.1 = read.csv('...\\Cohort1.csv',row.names = 1)
cohort.2 = read.csv('...\\Cohort2.csv',row.names = 1)

cohort.1$group <- factor(cohort.1$group, levels = c("B-enriched","T-enriched"))
cohort.2$group <- factor(cohort.2$group, levels = c("B-enriched","T-enriched"))

features <- setdiff(colnames(cohort.1), "group")
oob_tbl <- data.frame(Subset=character(), NumFeatures=integer(), OOB_AUC=numeric())
best_auc <- -Inf; best_sub <- NULL; best_model <- NULL; best_roc <- NULL

set.seed(222)
for (k in 1:length(features)) {
  subs <- combn(features, k, simplify = FALSE)
  for (sub in subs) {
    rf_fit <- randomForest(x = cohort.1[, sub, drop=FALSE],
                           y = cohort.1$group,
                           ntree = 500)
    oob_probs <- rf_fit$votes[, "T-enriched"]
    roc_oob   <- roc(cohort.1$group, oob_probs,
                     levels=c("B-enriched","T-enriched"), quiet=TRUE)
    auc_oob   <- as.numeric(auc(roc_oob))
    
    oob_tbl <- rbind(oob_tbl,
                     data.frame(Subset=paste(sub, collapse = "+"),
                                NumFeatures=length(sub),
                                OOB_AUC=auc_oob))
    if ( (auc_oob > best_auc + 1e-12) ||
         (abs(auc_oob - best_auc) <= 1e-12 && length(sub) < length(best_sub)) ) {
      best_auc  <- auc_oob
      best_sub  <- sub
      best_model<- rf_fit
      best_roc  <- roc_oob
    }
  }
}

oob_tbl <- oob_tbl[order(oob_tbl$NumFeatures), ]
oob_tbl$index = c(1:127)
write.csv(oob_tbl,'...\\oob tbl result.csv')
best_idx <- 42
a = ggplot(oob_tbl, aes(x = index, y = OOB_AUC)) +
  geom_point(color = "#2E86C1", size = 2) +
  geom_line(color = "#2E86C1", linewidth = 0.8) +
  geom_point(data = subset(oob_tbl, index == best_idx),
             aes(x = index, y = OOB_AUC),
             color = "red", size = 3) +
  theme_minimal(base_size = 14) +
  labs(x = "Feature Subset Index",
       y = "OOB AUC",
       title = "OOB AUC across Feature Subsets") +
  scale_x_continuous(breaks = seq(0, max(oob_tbl$index), by = 5))


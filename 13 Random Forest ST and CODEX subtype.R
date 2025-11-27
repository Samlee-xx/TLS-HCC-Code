
###### Random Forest ######
### Spatial Transcriptom
library(randomForest)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(scales)
set.seed(777)

setwd('...\\06 Cell2location')
data = read.csv('tls_cell2location_out.csv')
tls=readRDS('...\\tls.rds')
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

data.norm.tb = subset(data.norm,group %in% c('T-enriched',
                                             'B-enriched'))
colnames(data.norm.tb)
ggplot(data.norm.tb,aes(b.linkage,t.linkage,colour = group))+
  geom_point()

data.norm.tb$group <- as.factor(data.norm.tb$group)
features <- data.norm.tb[,c(1:7)]
label <- data.norm.tb$group

rf_model <- randomForest(x = features, y = as.factor(label), importance = TRUE, ntree = 500)
importance_df <- importance(rf_model)
importance_df <- data.frame(Feature = rownames(importance_df),
                            Importance = importance_df[, "MeanDecreaseGini"])
a = ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance (Random Forest)", x = "Feature", y = "Importance (Gini)") +
  theme_minimal()

pdf(file = '...\\ST featureImportance.pdf',
    width = 3,height = 3)
print(a)
dev.off()

features_cols <- colnames(features) 
means_by_group <- data.norm.tb %>%
  group_by(group) %>%
  summarise(across(all_of(features_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
stopifnot(all(c("B-enriched", "T-enriched") %in% means_by_group$group))

b_vec <- as.numeric(means_by_group[means_by_group$group == "B-enriched", features_cols])
t_vec <- as.numeric(means_by_group[means_by_group$group == "T-enriched", features_cols])
names(b_vec) <- features_cols
names(t_vec) <- features_cols

eps <- 1e-3
log2FC <- log2((b_vec + eps) / (t_vec + eps))

fc_df <- tibble(
  Feature   = features_cols,
  log2FC    = as.numeric(log2FC),
  enriched  = ifelse(log2FC > 0, "B-enriched", "T-enriched")
)

plot_df <- importance_df %>%
  inner_join(fc_df, by = "Feature") %>%
  arrange(Importance)
max_abs_fc <- max(abs(plot_df$log2FC), na.rm = TRUE)

a = ggplot(plot_df, aes(x = reorder(Feature, Importance), y = Importance, fill = log2FC)) +
  geom_col(color = "grey20", width = 0.8) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#43b0f1", mid = "white", high = "#e12729", midpoint = 0,
    limits = c(-max_abs_fc, max_abs_fc),  
    name = "log2 FC (B/T)",
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    title = "Feature Importance (RF) with Group Enrichment",
    x = "Feature", y = "Mean Decrease Gini"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )


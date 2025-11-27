
###### Cross K analysis of neutrophil and cholangiocyte ######
library(tidyverse)
library(spatstat)
library(tidyverse)
library(spatstat)
library(furrr)

tls.meta = readRDS('E:\\01-TLS\\09-CODEX Alignment\\07-HD TLS Quary\\Analysis\\03 Correlation\\TLS codex.rds')
colnames(tls.meta)
new_features <- tls.meta[, c("B", "Bile_Duct", "CD4T", "CD8T", "Endothelial", "Fibroblast", "Macrophage", "Neutrophil")]
new_features$location <- sub(".*_(N|P|T).*", "\\1", rownames(new_features))
table(new_features$group,new_features$location)
table(new_features$location)
ta.new.features = subset(new_features,location %in% c('P','T'))

rownames(ta.new.features)

codex.total = read.csv('E:\\01-TLS\\09-CODEX Alignment\\06-TLS DBSCAN\\03-Codex Metadata\\Total CODEX Metadata 0206.csv')
codex.total$tls.id = paste0(codex.total$slide,"_",codex.total$fov)
rownames(ta.new.features) %in% unique(codex.total$tls.id)
rownames(ta.new.features)[25:34] = c('239117_P_1',
                                     '239117_P_2',
                                     '239117_P_3',
                                     '239117_T_1',
                                     '239117_T_2',
                                     '239117_T_4',
                                     '239117_T_5',
                                     '239117_T_6',
                                     '239117_T_7',
                                     '239117_T_8')

codex.ta.tls = subset(codex.total,tls.id %in% rownames(ta.new.features))
codex.ta.tls$cellid = paste0(codex.ta.tls$fov,'_',codex.ta.tls$label)
table(codex.ta.tls$celltype)
ta.bile.duct <- codex.ta.tls[codex.ta.tls$celltype == 'Bile_Duct', 'cellid']
bile.cells <- codex.ta.tls[codex.ta.tls$cellid %in% ta.bile.duct, ]
codex.ta.tls$centroid.0 <- as.numeric(codex.ta.tls$centroid.0)
codex.ta.tls$centroid.1 <- as.numeric(codex.ta.tls$centroid.1)

target_r <- 100
df_cross_k <- codex.ta.tls %>%
  group_by(tls.id) %>%
  nest() %>%
  mutate(cross_k = future_map(data, function(dat) {
    bd <- dat %>% filter(celltype == "Bile_Duct") %>%
      select(centroid.0, centroid.1) %>% mutate(type = "Bile_Duct")
    neu <- dat %>% filter(celltype == "Neutrophil") %>%
      select(centroid.0, centroid.1) %>% mutate(type = "Neutrophil")
    if (nrow(bd) < 3 | nrow(neu) < 3) return(NULL)
    
    # merge data and create window
    combined_data <- bind_rows(bd, neu)
    combined_ppp <- ppp(
      x = combined_data$centroid.0,
      y = combined_data$centroid.1,
      window = owin(
        xrange = range(dat$centroid.0, na.rm = TRUE),
        yrange = range(dat$centroid.1, na.rm = TRUE)
      ),
      marks = as.factor(combined_data$type)
    )
    
    # calculate cross-type K function(Ripley isotropic correct)
    cross_k <- Kcross(
      combined_ppp,
      i = "Bile_Duct",
      j = "Neutrophil",
      correction = "isotropic",
      rmax = target_r
    )
    
    # find the closest target_r index
    idx <- which.min(abs(cross_k$r - target_r))
    
    # return tibble
    tibble(
      r = cross_k$r[idx],
      theo = cross_k$theo[idx],
      obs = cross_k$iso[idx],  # isotropic correction
      normalized_diff = (cross_k$iso[idx] - cross_k$theo[idx]) / cross_k$theo[idx]
    )
  })) %>%
  select(-data) %>%
  unnest(cols = cross_k) %>%
  drop_na(normalized_diff) 

a = ggplot(df_cross_k, aes(x = reorder(tls.id, normalized_diff), y = normalized_diff)) +
  geom_segment(aes(xend = tls.id, yend = 0), color = "gray80", linetype = "dotted") +
  geom_point(aes(color = normalized_diff > 0), size = 1.8) +
  geom_hline(yintercept = 0, color = "black") +
  scale_color_manual(values = c("TRUE" = "#adc178", "FALSE" = "#ffa0c5")) +
  labs(
    title = "Cross-K (r = 100 px, ???25 µm): Bile Duct vs Neutrophil",
    x = "TLS ID",
    y = "Normalized Difference (Observed - Theoretical) / Theoretical"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

wilcox_res <- wilcox.test(df_cross_k$normalized_diff, mu = 0)
print(wilcox_res)



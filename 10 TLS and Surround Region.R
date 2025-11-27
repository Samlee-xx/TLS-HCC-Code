
###### TLS and TLS surround ######
### 01 Paired comparision ----
# total cellmetadata
total.cell = readRDS('...\\01 Total Cell Metadata-0318.rds')
table(total.cell$category) 
colnames(total.cell)
ta.tls.meta = readRDS('...\\TA_TLS_Metadata.rds')
ta.tls.fov = unique(ta.tls.meta$fov)

ta.tls.fov.meta = subset(total.cell,fov %in% ta.tls.fov)
table(ta.tls.fov.meta$new.celltype)
ta.tls.fov.meta = subset(ta.tls.fov.meta, !(new.celltype %in% "Tumor"))

library(data.table)
library(RANN)

setDT(ta.tls.fov.meta)
ta.tls.fov.meta <- ta.tls.fov.meta[!is.na(centroid.0) & !is.na(centroid.1)]
ta.tls.fov.meta[, is_tls := !is.na(tls_id) & tls_id != "Out TLS Cells"]

classify_one_fov <- function(DT, cutoff = 100) {
  DT <- data.table::copy(DT) 
  
  if (!any(DT$is_tls)) {
    DT[, `:=`(nn_dist_to_tls = NA_real_,
              region = fifelse(is_tls, "TLS", "other cells"))]
    return(DT)
  }
  
  tls_xy <- as.matrix(DT[is_tls == TRUE, .(centroid.0, centroid.1)])
  
  q_idx <- which(DT$is_tls == FALSE)
  if (length(q_idx)) {
    nn <- nn2(
      data  = tls_xy,
      query = as.matrix(DT[q_idx, .(centroid.0, centroid.1)]),
      k = 1
    )
    DT[q_idx, nn_dist_to_tls := nn$nn.dists[, 1]]
  }
  
  DT[is_tls == TRUE, nn_dist_to_tls := 0]
  
  DT[, region := fifelse(
    is_tls,
    "TLS",
    fifelse(!is.na(nn_dist_to_tls) & nn_dist_to_tls <= cutoff,
            "surrounding TLS", "other cells")
  )]
  
  DT
}

res <- ta.tls.fov.meta[, classify_one_fov(.SD, cutoff = 200), by = fov]
table(res$region)

res = subset(res, !(new.celltype %in% "Unassign"))
surround.tls = subset(res,region == 'surrounding TLS')
used.tls.meta = subset(res,region == 'TLS')

back.immune.pro = table(surround.tls$new.celltype) %>% as.data.frame()
tls.immune.pro = table(used.tls.meta$new.celltype) %>% as.data.frame()
back.immune.pro$Freq2 = tls.immune.pro$Freq
colnames(back.immune.pro) = c('celltype','Backgroud','TLS')
rownames(back.immune.pro) = back.immune.pro$celltype

back.immune.pro = melt(back.immune.pro)

colnames(back.immune.pro) = c('Celltype','Location','Freq')
library(plyr)
library(ggplot2)
library(ggplot2)
library(ggprism)
library(reshape)
library(ggalluvial)
library(RColorBrewer)
df.nt=ddply(back.immune.pro,'Location',transform,
            percentage = Freq/sum(Freq))

colorlist <- c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f",
               "#b2db87","#7ee7bb","#64cccf","#a9dce6","#a48cbe",'#e4b7d6','#757575','gray')

a= ggplot(df.nt, aes( x = Location,y=percentage,fill = Celltype,
                      stratum = Celltype, alluvium = Celltype))+
  geom_stratum(width = 0.5, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear")+
  scale_fill_manual(values =colorlist)+
  theme_classic()


# each celltype 
# keep TLS & surrounding TLS region
res2 <- res[region %in% c("TLS", "surrounding TLS")]

fovs  <- unique(res2$fov)
ctypes <- unique(res2$new.celltype)
regions <- c("TLS", "surrounding TLS")

# fov × region × celltype
cnt <- res2[, .N, by = .(fov, region, new.celltype)]

cnt_full <- CJ(fov = fovs, region = regions, new.celltype = ctypes)
cnt_full <- cnt_full[cnt, on = .(fov, region, new.celltype)]
cnt_full[is.na(N), N := 0L]

# FOV × region 
denom <- cnt_full[, .(region_total = sum(N)), by = .(fov, region)]

prop_dt <- cnt_full[denom, on = .(fov, region)]
prop_dt[, percentage := ifelse(region_total > 0, 100 * N / region_total, 0)]

final.df <- prop_dt[, .(fov, celltype = new.celltype,
                        variable = factor(region, levels = c("TLS","surrounding TLS")),
                        value = percentage)]

keep_fov <- final.df[, .N, by = .(fov, variable)][, .N, by = fov][N == 2, fov]
final.df <- final.df[fov %in% keep_fov]

Name <- "Surrounding vs TLS (T-B)"
cell_of_interest <- "T-B" # cell types 

table(final.df$celltype)
plot_df <- final.df[celltype == cell_of_interest]

keep_fov <- plot_df[, uniqueN(variable), by = fov][V1 == 2, fov]
plot_df <- plot_df[fov %in% keep_fov]

a = ggpaired(plot_df, x = "variable", y = "value", fill = "variable",
             id = "fov",point.size = 1, line.alpha = 0.01,
             add = "jitter", line.color = "#e2e2e4", line.size = 0.1,point.alpha =0.1,
             palette = c("#43b0f1", "#d74a49"),
             xlab = " ", ylab = "Percentage", title = Name,
             legend.title = " ", show.legend = FALSE) +
  stat_compare_means(method = "wilcox.test",
                     paired = TRUE,
                     label = "p.format") +    
  theme(legend.position = "none")

### 02 spatial enrichment ----
total.cell = readRDS('...\\01 Total Cell Metadata-0318.rds')
table(total.cell$category) 
colnames(total.cell)
ta.tls.meta = readRDS('...\\TA_TLS_Metadata.rds')
ta.tls.fov = unique(ta.tls.meta$fov) 

ta.tls.fov.meta = subset(total.cell,fov %in% ta.tls.fov)
table(ta.tls.fov.meta$new.celltype)
ta.tls.fov.meta = subset(ta.tls.fov.meta, !(new.celltype %in% "Tumor"))

library(data.table)
library(RANN)

setDT(ta.tls.fov.meta)
ta.tls.fov.meta <- ta.tls.fov.meta[!is.na(centroid.0) & !is.na(centroid.1)]
ta.tls.fov.meta[, is_tls := !is.na(tls_id) & tls_id != "Out TLS Cells"]

classify_one_fov <- function(DT, cutoff = 100) {
  DT <- data.table::copy(DT) 
  
  if (!any(DT$is_tls)) {
    DT[, `:=`(nn_dist_to_tls = NA_real_,
              region = fifelse(is_tls, "TLS", "other cells"))]
    return(DT)
  }
  
  tls_xy <- as.matrix(DT[is_tls == TRUE, .(centroid.0, centroid.1)])
  
  q_idx <- which(DT$is_tls == FALSE)
  if (length(q_idx)) {
    nn <- nn2(
      data  = tls_xy,
      query = as.matrix(DT[q_idx, .(centroid.0, centroid.1)]),
      k = 1
    )
    DT[q_idx, nn_dist_to_tls := nn$nn.dists[, 1]]
  }
  
  DT[is_tls == TRUE, nn_dist_to_tls := 0]
  
  DT[, region := fifelse(
    is_tls,
    "TLS",
    fifelse(!is.na(nn_dist_to_tls) & nn_dist_to_tls <= cutoff,
            "surrounding TLS", "other cells")
  )]
  
  DT
}

res <- ta.tls.fov.meta[, classify_one_fov(.SD, cutoff = 200), by = fov]

table(res$region)

res = subset(res,region %in% c('TLS','surrounding TLS'))

celltype.to.analyze <- "Macrophage"
table(res$new.celltype)
neutro_data <- res %>%
  filter(new.celltype == celltype.to.analyze)

neutro_fov_count <- neutro_data %>%
  dplyr::group_by(fov) %>%
  dplyr::summarise(n_celltype = dplyr::n()) %>%
  dplyr::filter(n_celltype > 50)

filtered_data <- neutro_data %>%
  filter(fov %in% neutro_fov_count$fov)

library(dplyr)
neutro_tls_stats <- filtered_data %>%
  dplyr::mutate(location = ifelse(region == "surrounding TLS", "Out", "In")) %>%
  dplyr::group_by(fov, location) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = location, values_from = n, values_fill = 0) %>%
  mutate(
    total = In + Out,
    in_ratio = In / total,
    out_ratio = Out / total
  )

tls_area_df <- res %>%
  dplyr::group_by(fov) %>%
  dplyr::summarise(
    total_cells = n(),
    tls_cells = sum(tls_id != "Out TLS Cells", na.rm = TRUE) 
  ) %>%
  dplyr:: mutate(tls_area_ratio = tls_cells / total_cells)

library(dplyr)
neutro_tls_stats <- neutro_tls_stats %>%
  left_join(dplyr::select(tls_area_df, fov, tls_area_ratio), by = "fov") %>%
  mutate(enrichment = in_ratio / tls_area_ratio)

a = ggplot(neutro_tls_stats, aes(x = "", y = enrichment)) +
  geom_violin(fill = '#8162e9', alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  theme_classic(base_size = 14) +
  labs(
    title = paste("Enrichment of", celltype.to.analyze, "in TLS"),
    subtitle = paste0("Wilcoxon test p = ", signif(wilcox.test(neutro_tls_stats$enrichment, mu = 1)$p.value, 3)),
    y = "Enrichment Score (Observed / Expected)",
    x = ""
  )+ylim(0.5,2)




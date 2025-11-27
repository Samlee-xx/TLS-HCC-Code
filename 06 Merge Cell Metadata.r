
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
cd4.t = readRDS('...\\01 CD4T\\01 CD4T.RDS')
cd8.t = readRDS('...\\02 CD8T\\01 CD8T.RDS')
neu = readRDS('...\\03 Neu\\01 Neu.RDS')
macro = readRDS('...\\04 Macro\\01 Macro.rds')
b = readRDS('...\\05 B - New verison\\01 B cell.rds')
stroma = readRDS('...\\01 Stroma 0809.rds')
stroma = readRDS('...\\06 Stroma\\01 Stroma.rds')
table(stroma$celltype)

total.meta = readRDS('E...\\01 Total Cell Metadata-0603.rds')
table(total.meta$celltype)
sum(is.na(total.meta$celltype)) 

cd4.t$celltype = as.character(cd4.t$celltype)
total.meta$subtype[match(colnames(cd4.t),total.meta$cellid)] = cd4.t$celltype 
cd8.t$celltype = as.character(cd8.t$celltype)
total.meta$subtype[match(colnames(cd8.t),total.meta$cellid)] = cd8.t$celltype
neu$celltype = as.character(neu$celltype)
total.meta$subtype[match(colnames(neu),total.meta$cellid)] = neu$celltype
macro$celltype = as.character(macro$celltype)
total.meta$subtype[match(colnames(macro),total.meta$cellid)] = macro$celltype
b$celltype = as.character(b$celltype)
total.meta$subtype[match(colnames(b),total.meta$cellid)] = b$celltype
stroma$celltype = as.character(stroma$celltype)
total.meta$subtype[match(colnames(stroma),total.meta$cellid)] = stroma$celltype
table(total.meta$subtype)
total.meta$subtype
total.meta$subtype = ifelse(total.meta$celltype == 'DC','DC',total.meta$subtype)
total.meta$subtype = ifelse(total.meta$celltype == 'T-B','T-B',total.meta$subtype)
table(total.meta$subtype)
table(total.meta$celltype)
total.meta$celltype = ifelse(total.meta$celltype == 'Artifact','Unassign',total.meta$celltype)
total.meta$subtype = ifelse(total.meta$subtype == 'Artifact','Unassign',total.meta$subtype)
total.meta$celltype[is.na(total.meta$celltype)] ='Unassign'
total.meta$subtype[is.na(total.meta$subtype)] ='Unassign'
table(total.meta$subtype)
table(total.meta$celltype)

setwd('...\\03 Merge')
saveRDS(total.meta,'01 Total Cell Metadata-0811.rds')


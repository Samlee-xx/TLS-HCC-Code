
###### Pre-process HD TLS Data ######

library(Seurat)
library(dplyr)
library(reshape2)
library(dplyr)
library(IOBR)
library(sva)
library(stringr)

hd_135121 = read.csv('...\\135121 TLS Barcode.csv')
barcode_135121 = Read10X('...\\square_008um\\filtered_feature_bc_matrix')

hd_239117 = read.csv('...\\239117 TLS Barcode.csv')
barcode_239117 = Read10X('...\\square_008um\\filtered_feature_bc_matrix')

hd_238966 = read.csv('...\\238966 TLS Barcode.csv')
barcode_238966 = Read10X('...\\square_008um\\filtered_feature_bc_matrix')

hd_178031 = read.csv('...\\178031 TLS Barcode.csv')
barcode_178031 = Read10X('...\\square_008um\\filtered_feature_bc_matrix')

hd_164697 = read.csv('...\\164697 TLS Barcode.csv')
barcode_164697 = Read10X('...\\square_008um\\filtered_feature_bc_matrix')

# one
hd.135121.filter = barcode_135121[,hd_135121$Barcode]
colnames(hd.135121.filter) == hd_135121$Barcode
hd.135121.filter = hd.135121.filter %>%as.data.frame()%>%t() %>% as.data.frame() 
hd.135121.filter$tls.id = hd_135121$TLS
hd.135121.expr = hd.135121.filter %>% group_by(tls.id) %>% summarise_each(funs = sum)
hd.135121.expr = hd.135121.expr%>% as.data.frame()
rownames(hd.135121.expr) = hd.135121.expr$tls.id
hd.135121.expr = hd.135121.expr[,-1]
hd.135121.expr = hd.135121.expr %>%t() %>% as.data.frame()
# two
hd.239117.filter = barcode_239117[,hd_239117$Barcode]
colnames(hd.239117.filter) == hd_239117$Barcode
hd.239117.filter = hd.239117.filter %>% as.data.frame() %>% t() %>% as.data.frame() 
hd.239117.filter$tls.id = hd_239117$TLS
hd.239117.expr = hd.239117.filter %>% group_by(tls.id) %>% summarise_each(funs = sum)
hd.239117.expr = hd.239117.expr%>% as.data.frame()
rownames(hd.239117.expr) = hd.239117.expr$tls.id
hd.239117.expr = hd.239117.expr[,-1]
hd.239117.expr = hd.239117.expr %>%t() %>% as.data.frame()
#three
hd.238966.filter = barcode_238966[,hd_238966$Barcode]
colnames(hd.238966.filter) == hd_238966$Barcode
hd.238966.filter = hd.238966.filter %>% as.data.frame() %>% t() %>% as.data.frame()
hd.238966.filter$tls.id = hd_238966$TLS
hd.238966.expr = hd.238966.filter %>% group_by(tls.id) %>% summarise_each(funs = sum)
hd.238966.expr = hd.238966.expr%>% as.data.frame()
rownames(hd.238966.expr) = hd.238966.expr$tls.id
hd.238966.expr = hd.238966.expr[,-1]
hd.238966.expr = hd.238966.expr %>%t() %>% as.data.frame()
#four 
hd.178031.filter = barcode_178031[,hd_178031$Barcode]
colnames(hd.178031.filter) == hd_178031$Barcode
hd.178031.filter = hd.178031.filter %>% as.data.frame() %>% t() %>% as.data.frame() 
hd.178031.filter$tls.id = hd_178031$TLS
hd.178031.expr = hd.178031.filter %>% group_by(tls.id) %>% summarise_each(funs = sum)
hd.178031.expr = hd.178031.expr%>% as.data.frame()
rownames(hd.178031.expr) = hd.178031.expr$tls.id
hd.178031.expr = hd.178031.expr[,-1]
hd.178031.expr = hd.178031.expr %>%t() %>% as.data.frame()
#five
hd.164697.filter = barcode_164697[,hd_164697$Barcode]
colnames(hd.164697.filter) == hd_164697$Barcode
hd.164697.filter = hd.164697.filter %>% as.data.frame() %>% t() %>% as.data.frame() 
hd.164697.filter$tls.id = hd_164697$TLS
hd.164697.expr = hd.164697.filter %>% group_by(tls.id) %>% summarise_each(funs = sum)
hd.164697.expr = hd.164697.expr%>% as.data.frame()
rownames(hd.164697.expr) = hd.164697.expr$tls.id
hd.164697.expr = hd.164697.expr[,-1]
hd.164697.expr = hd.164697.expr %>%t() %>% as.data.frame()

gene_length.fliter = read.csv('E:\\01-TLS\\09-CODEX Alignment\\07-HD TLS Quary\\All_hg19gene_len.csv')
TPM.135121=count2tpm(as.matrix(hd.135121.expr),effLength = gene_length.fliter,
                     id='Gene',
                     length = 'Length',
                     gene_symbol = 'Gene')

TPM.239117=count2tpm(as.matrix(hd.239117.expr),effLength = gene_length.fliter,
                     id='Gene',
                     length = 'Length',
                     gene_symbol = 'Gene')

TPM.238966=count2tpm(as.matrix(hd.238966.expr),effLength = gene_length.fliter,
                     id='Gene',
                     length = 'Length',
                     gene_symbol = 'Gene')

TPM.178031=count2tpm(as.matrix(hd.178031.expr),effLength = gene_length.fliter,
                     id='Gene',
                     length = 'Length',
                     gene_symbol = 'Gene')

TPM.164697=count2tpm(as.matrix(hd.164697.expr),effLength = gene_length.fliter,
                     id='Gene',
                     length = 'Length',
                     gene_symbol = 'Gene')

save(list = c('TPM.135121','TPM.164697','TPM.178031','TPM.238966','TPM.239117'),
     file = '...\\HD TPM.rdata')



setwd('...\\03 Keep used TLS')
total.meta = readRDS('...\\01 Total Cell Metadata-0603.rds')

library(stringr)
back.info = str_split(total.meta$cellid,'_',simplify = TRUE)
back.info = back.info %>% as.data.frame()
back.info$V4 = ifelse(back.info$V4 == 'B','Backgroud','TLS')
total.meta$location = back.info$V2
table(total.meta$location)

table(total.meta$slide)
total.meta$stage = case_when(total.meta$slide =='Sample_123942'~'Early',
                             total.meta$slide =='Sample_124989'~'Early',
                             total.meta$slide =='Sample_133698'~'Late',
                             total.meta$slide =='Sample_134024'~'Middle',
                             total.meta$slide =='Sample_134205'~'Early',
                             total.meta$slide =='Sample_136254'~'Middle',
                             total.meta$slide =='Sample_136759'~'Middle',
                             total.meta$slide =='Sample_176552'~'Late',
                             total.meta$slide =='Sample_86130'~'Early')

table(total.meta$subtype)
total.meta = subset(total.meta,Backgroud =='TLS') 
total.meta = filter(total.meta,cluster_label >-1) 

remove.fov = c('123942_N_10','123942_N_21','123942_N_30','123942_T_13','123942_T_16','123942_T_19',
               '123942_T_22',
               '124989_I_2','124989_N_1','124989_N_2','124989_N_10','124989_N_13','124989_T_5',
               '124989_T_8','124989_T_9','124989_T_10',
               '133698_T_2','133698_T_3','133698_T_4',
               '134024_T_9','134024_T_10',
               '134205_T_2',
               '136759_T_1',
               '176552_L_4')

total.meta = subset(total.meta, !(fov %in% remove.fov))
colnames(total.meta)

subtype.table = total.meta[,c('tls_id','celltype')]
subtype.table = dcast(subtype.table,tls_id~celltype)
rownames(subtype.table) = subtype.table$tls_id
colnames(subtype.table)
subtype.table = subtype.table[,-c(1)]

normalize = function(x){
  x/sum(x)
}

subtype.table = apply(subtype.table, 1, normalize)
subtype.table = subtype.table %>% t() %>% as.data.frame()
subtype.table = subtype.table*100

subtype.table = filter(subtype.table,Unassign <=30)
subtype.table = filter(subtype.table,B+CD4T+CD8T+`T-B` >=10)
subtype.table = filter(subtype.table,Stroma <=60)
keep.tls = rownames(subtype.table)
saveRDS(keep.tls,'01 Used TLS name-V0603.rds')





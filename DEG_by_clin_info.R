
# root - grp_dir - files ...
# 구조의 데이터 RNAseq

setwd("C:/Users/user/Desktop/git_manager/RNA-seq_ej")

library(dplyr)

getwd()

# root_dir = 'E:/stemcell/RNAseq/gdc/count/origin'
# root_dir = 'E:/stemcell/RNAseq/gdc/count/mt_mut'
# root_dir = 'E:/stemcell/RNAseq/gdc/count/Diff_whole'
root_dir = 'E:/stemcell/RNAseq/gdc/count/D_Age'
# root_dir = 'E:/stemcell/RNAseq/gdc/count/Origin_comp/AM_Blo'
# root_dir = 'E:/stemcell/RNAseq/gdc/count/Origin_comp/Blo_Fib'

# save_path = 'table_data/mt_mut_table.tsv'

# output_f_name = 'mtMut_DEG'
# output_f_name = 'Diff_whole_DEG'
output_f_name = 'AGE_DEG'
# output_f_name = 'AM_Blo_DEG'
# output_f_name = 'AM_Fib_DEG'
# output_f_name = 'Blo_Fib_DEG'

# target_sample_grp_name = 'MtMut'
# base_sample_grp_name = 'WildType'
target_sample_grp_name = 'Old'
base_sample_grp_name = 'Fetal'
# target_sample_grp_name = 'Blo'
# base_sample_grp_name = 'AM'
# target_sample_grp_name = 'Fib'
# base_sample_grp_name = 'AM'
# target_sample_grp_name = 'Fib'
# base_sample_grp_name = 'Blo'
# save_path = 'table_data/origin_mut_table.tsv'

# root_dir = 'C:/AMC_proj/stemcell/RNAseq/RNAseq/gdc/count/mt_mut'
# save_path = 'table_data/mt_mut_table.tsv'

old.path <- setwd(root_dir)
getwd()


count_file_lst <- list.files(recursive = T, pattern = '*.count')
count_file_lst
count_file_lst[1]



def_make_DEG_inputs <- function(file_lst){
  
  
  first.sample.n = read.delim(file_lst[1],
                              header = F, row.names = 1)
  
  # head(first.sample.n)
  # class(first.sample.n) # df
  # 
  count.table = data.frame(first.sample.n)
  # head(count.table)
  # str(count.table)
  # class(count.table)
  # 
  # 
  # # gsub(찾을 것, 바꿀 것, 열 지정)
  # file_lst
  # # name_list = gsub(".count","", file_lst)
  # # name_list
  
  grp_vec <- c()
  sam_name_vec <- c()
  
  for (f_tag in file_lst){
    split_vec1 <- strsplit(f_tag, split = '/')
    grp_ <- split_vec1[[1]][1]
    f_name <- split_vec1[[1]][2]
    f_name = gsub(".count","", f_name)
    
    # print(grp_)
    # print(f_name)
    
    grp_vec <- c(grp_vec, grp_)
    sam_name_vec <- c(sam_name_vec, f_name)
  }
  
  print(grp_vec)
  print(sam_name_vec)
  
  # library(tibble)
  
  for (i in file_lst[2:length(file_lst)]) {
    column.n = read.delim(i, header = F, row.names = 1)
    # print(tibble(column.n))
    count.table = cbind(count.table, s = column.n)
  }
  
  colnames(count.table) = sam_name_vec
  
  head(count.table)
  
  
  return_lst <- list(count.table,
                     grp_vec,
                     sam_name_vec)
  
  names(return_lst) <- c('count_table', 'grp_vec', 'sam_name_vec')
  
  return(return_lst)
}

input_lst <- def_make_DEG_inputs(count_file_lst)
input_lst

count.table <- input_lst$count_table
grp.vec <- input_lst$grp_vec
grp.vec <- as.factor(grp.vec)
sam.name.vec <- input_lst$sam_name_vec

count.table
grp.vec
sam.name.vec
# write.csv(count.table, file=file.path(save_dir, save_csv_name))
# write.table(count.table, file=file.path(save_dir, save_tsv_name), sep = "\t")
# write.table(count.table, file=file.path(save_dir, save_gsea_name), sep = "\t")
# write.table(count.table, file=save_path, sep = "\t")


# DEG 분석


library(edgeR)

library(tibble)


def_DEG_analysis_glm <- function(count_table_, grp_vec_, sam_name_vec_,
                                 target_sample_count, base_sample_count,
                                 target_samplegrp_name, base_samplegrp_name){
  
  meta.data <- tibble(data.frame(sam_name_vec_, grp_vec_))
  names(meta.data) <- c('SampleID', 'Condition')
  meta.data$Condition <- factor(meta.data$Condition)
  
  print(meta.data, n=dim(meta.data)[1])
  print(str(meta.data))
  print(levels(meta.data$Condition))
  
  
  
  design <- model.matrix(~Condition, data=meta.data)
  print(design)

  y <- DGEList(counts = count_table_,
               genes = rownames(count_table_),
               group = grp_vec_)
  
  
  
  keep <- filterByExpr(y, group = grp_vec_)
  ?filterByExpr
  print(table(keep))
  
  y <- y[keep, , keep.lib.sizes=FALSE]

  y <- calcNormFactors(y)
  y <- estimateGLMRobustDisp(y, design = design)

  print(y$samples)
  #
  # fit <- glmFit(y, design)

  MDS_data <- plotMDS(y, top=20000)
  
  plotMDS(y, top=20000, labels=c(rep(target_samplegrp_name, target_sample_count), 
                                 rep(base_samplegrp_name, base_sample_count)),
          col = c(rep("tomato1",target_sample_count), 
                  rep("steelblue1",base_sample_count)))
  # plot_data <-

  fit <- glmFit(y, design = design)
  lrt <- glmLRT(fit, coef = 2)
  res.table <- topTags(lrt, n=dim(count_table_)[1], sort.by = 'none')$table
  
  plotMD(lrt)
  abline(h=c(-1.5, 1.5), col='blue')

  return_lst <- list(y, lrt, design, res.table, count_table_)
  names(return_lst) <- c('y', 'lrt', 'design', 'result_table', 
                         'raw_count_table')

  return(return_lst)
  
}

# ?glmFit
# ?glmLRT

base_count <- sum(grp.vec == levels(grp.vec)[1])
target_count <- sum(grp.vec == levels(grp.vec)[2])
target_count
base_count
deg.result.list <- def_DEG_analysis_glm(count.table, grp.vec, sam.name.vec,
                                        target_count, base_count,
                                        target_sample_grp_name,
                                        base_sample_grp_name)

deg.result.list$lrt
deg.result.list$y
deg.result.list$y$AveLogCPM
summary(deg.result.list$raw_count_table[,1])
summary(deg.result.list$y$counts[,1])


# test.cpm <- cpm(deg.result.list$y, normalized.lib.sizes = T, log = T)
# test.cpm


nrow(deg.result.list$raw_count_table)
deg.result.list

sum(deg.result.list$result_table$PValue<0.01) # 303
sum(deg.result.list$result_table$PValue<0.05) # 1128

sum(deg.result.list$result_table$FDR<0.01) # 43
sum(deg.result.list$result_table$FDR<0.05) # 50

# 
# # TMM값 붙여주기  / 나중에 꼭 필요하면 코드 수정 (cpm 커트해서 매치가 안됨됨)
# ?cpm
# ?rpkm
# 
# TMM <- cpm(deg.result.list$y, normalized.lib.sizes = T, log = T)
# head(TMM)
# TMM_colName <- paste("TMMValue_", colnames(TMM), sep = '')
# TMM_colName
# colnames(TMM) <- TMM_colName
# TMM
# 
# TMM.mean <- apply(TMM, MARGIN = 1, FUN = mean)
# class(TMM)
# TMM <- as.data.frame(TMM)
# TMM$TMM.mean <- TMM.mean
# 
# # raw count 값 붙여주기
# raw_cnt_data <- deg.result.list$raw_count_table
# raw_cnt_data
# 
# RawCount_colName <- paste("RawCount_", 
#                           colnames(deg.result.list$raw_count_table), sep = '')
# colnames(raw_cnt_data) <- RawCount_colName
# raw_cnt_data
# 
# deg.result.list$result_table
# 
# save_DEG_result <- cbind(deg.result.list$result_table,
#                          TMM,
#                          raw_cnt_data)
# 
# library(writexl)
# 
# write_xlsx(save_DEG_result, path = paste0(output_f_name, '_result.xlsx'))
# 
# 
# getwd()



deg.result.list$result_table[deg.result.list$result_table$FDR<0.01, ]

deg.result.list$result_table[deg.result.list$result_table$FDR<0.05, ]





write.table(deg.result.list$result_table, file=paste0(output_f_name, '_for_Vplot.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(deg.result.list$result_table, file=paste0(output_f_name, '_for_Vplot.csv'), quote = F, row.names = F)




res.filtered <- subset(deg.result.list$result_table, 
                           abs(logFC)>=1.5 & PValue<0.05)
dim(res.filtered) 


print('Filtering Result')
head(res.filtered)
dim(res.filtered)


up_regulated <- subset(res.filtered, logFC > 0)
down_regulated <- subset(res.filtered, logFC < 0)

# up_regulated
# down_regulated

up_genes <- up_regulated$genes
down_genes <- down_regulated$genes

library(writexl)

print('up/down reg gene count')
length(up_genes)
length(down_genes)


write.table(res.filtered, file=paste0(output_f_name, '_filtered.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(res.filtered, file=paste0(output_f_name, '_filtered.csv'), quote = F, row.names = F)




write.table(up_regulated, file=paste0(output_f_name, '_upreg.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(up_regulated, file=paste0(output_f_name, '_upreg.csv'), quote = F, row.names = F)

write.table(up_genes, file=paste0(output_f_name, '_upreg_genes.txt'), sep='\n',quote = FALSE,
            row.names = F, col.names = F)

write_xlsx(up_regulated, path = paste0(output_f_name, '_upreg.xlsx'))



write.table(down_regulated, file=paste0(output_f_name, '_downreg.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(down_regulated, file=paste0(output_f_name, '_downreg.csv'), quote = F, row.names = F)

write.table(down_genes, file=paste0(output_f_name, '_downreg_genes.txt'), sep='\n',quote = FALSE,
            row.names = F, col.names = F)

write_xlsx(down_regulated, path = paste0(output_f_name, '_downreg.xlsx'))






# 
# # for heatmap
# 
# class(deg.result.list$raw_count_table)
# deg.result.list$y$counts
# DEG.raw.table <- deg.result.list$raw_count_table
# DEG.raw.table
# DEG.raw.table$genes <- rownames(DEG.raw.table)
# 
# DEG.raw.table <- DEG.raw.table %>% relocate(genes)
# DEG.raw.table
# 
# 
# 
# filtered.genes <- res.filtered %>% select(genes)
# class(filtered.genes)
# filtered.genes
# 
# filtered_data_for_heatmap <- inner_join(filtered.genes, DEG.raw.table, by='genes')
# filtered_data_for_heatmap
# 
# library(tidyverse)
# 
# to_mat <- filtered_data_for_heatmap %>% 
#   remove_rownames %>% 
#   column_to_rownames(var="genes") %>% 
#   as.matrix()
# 
# to_mat
# class(to_mat)
# # to_mat[to_mat['29-B-Beta-cell']>0,]
# # sum(to_mat[,'29-B-Beta-cell']>=0)
# dim(to_mat)



# TMM <- cpm(to_mat, normalized.lib.sizes = T, log = T)
# head(TMM)
# TMM
# 
# tmp_y <- deg.result.list$y
# tmp_y
# tmp_y_cpm <- cpm(tmp_y, normalized.lib.sizes = T, log = T)
# class(tmp_y_cpm)
# sum(tmp_y_cpm[, '29-B-Beta-cell'] >10)
# 
# class(tmp_y$genes) # df
# class(tmp_y$counts) # mat
# 
# tmp_y$counts <- to_mat
# tmp_y
# tmp_y$genes <- filtered.genes
# tmp_y



# 
# deg.result.list.filterd <- def_DEG_analysis_glm(to_mat, grp.vec, sam.name.vec,
#                                         target_count, base_count,
#                                         target_sample_grp_name,
#                                         base_sample_grp_name)

# ?cpm
# TMM.heat <- cpm(deg.result.list.filterd$y, normalized.lib.sizes = T, log = T)




# 
# TMM.heat <- cpm(to_mat, normalized.lib.sizes = T, log = T)
# TMM.heat
# # 
# # deg.result.list.filterd$y
# # TMM.heat
# 
# 
# library(pheatmap)
# 
# pheatmap(TMM.heat, cluster_rows = T , fontsize_row = 3.5,
#          cluster_cols = F, scale = 'row')
# 
# pheatmap(TMM.heat, cluster_rows = T , fontsize_row = 3.5,
#          cluster_cols = F)












































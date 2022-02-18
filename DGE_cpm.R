

library(edgeR)

data_dir = "E:/stemcell/RNAseq/gdc/count/mt_mut/clone_comp/Tera/table_data"
table_name = "Tera_count.tsv"

output_f_name = 'Tera_DEG_result'

group_vec <- c(1, 2)

no_replicate_bcv <- 0.4

logFC_filter <- 1.5 # abs(logFC)
p_value <- 0.05



old.path <- setwd(data_dir)
getwd()


counts <- read.delim(file.path(data_dir, table_name), row.names = 1)
# head(counts)
# dim(counts)

raw_count <- rowSums(counts) > 0 # 최소 한쪽에서 발현이 0보다 많이 된것
# raw_count
raw.over.0 <- counts[raw_count,]
print('raw count (at least one)')
dim(raw.over.0)

samples <-  colnames(counts)
# samples
# counts[counts$hiPS29.A.EB==0,]
# sum(counts$X29.A.Beta.cell) # 118843204
# sum(counts$X29.B.Beta.cell) # 98828820
?cpm
cpm_count <- cpm(counts)
# cpm_count

# 트리밍할 유전자 정해주는 과정
thresh <- cpm_count > 0.3 # (y = (총 합 * 0.5) / 백만)
                          
# head(thresh)

print('CPM 기준: sum(bool)')
table(rowSums(thresh))

keep <- rowSums(thresh) >= 2 # 둘 다 역치 넘는것만 킵
# keep <- rowSums(thresh) >= 1 # 적어도 한개일때 킵
print('True = 둘 다 기준 이상')
summary(keep)
# head(keep)

# head(counts)

counts.keep <- counts[keep,] # keep에 해당하는 행(gene)만 뽑기
# head(counts.keep)

# counts.keep[counts.keep$X29.A.P50.Tera<10, ]
# counts.keep[counts.keep$X29.B.P50.Tera<10, ]

# dim(counts.keep) # 숫자 일치 확인

# sum(counts.keep$X29.A.P50.Tera) # 118843204
# sum(counts.keep$X29.B.P50.Tera) # 98828820
# counts.keep


d0 <- DGEList(counts.keep, group = group_vec) # colSums(counts.keep)이 lib.size임
# count_df와 sample_df로 이루어 져 있다.
# head(d0)
# class(d0)
# str(d0)

# head(d0$counts[, 'hiPS29.A.EB'])
d0$samples


# norm.factors 계산 및 추가
?calcNormFactors
y = calcNormFactors(d0) 
head(y)
# y$samples

# plotMDS(y)

cpm_count <- cpm(y, log=T) # log CPM 값
# head(cpm_count)

# class(cpm_count)



# 분산 추정 과정 
## no rep일때는 명시적으로 정해줘야함

# y1 = estimateDisp(y) #  rep sample이 없을때 : (경고) There is no replication

# y




#  exactTest
# ?exactTest

# no rep일때
# ?exactTest
et_no_rep <- exactTest(y, dispersion = no_replicate_bcv^2)
# et_no_rep
# class(et_no_rep)

et_no_rep$genes <- rownames(et_no_rep$table)
# et_no_rep$genes

# dim(et_no_rep)
# str(et_no_rep)

# dim(et_no_rep$table)
row_len <- dim(et_no_rep$table)[1]
# row_len

top50 <- topTags(et_no_rep, n = 50)
top50_table <- top50$table

# plotMD(top50)


top_res <- topTags(et_no_rep, n = row_len) # fdr 값 새겨진 table get
top_res_table <- top_res$table
# dim(top_res_table)

library(tibble)

top_res_table <- rownames_to_column(top_res_table, 'Gene')
top50_res_table <- rownames_to_column(top50_table, 'Gene')

head(top_res_table)
dim(top_res_table)

# head(top50_res_table)
# dim(top50_res_table)



# top_res_filtered <- subset(top_res_table, 
#                              abs(logFC)>logFC_filter & PValue<=p_value &
#                                FDR<=0.05)
top_res_filtered <- subset(top_res_table, 
                           abs(logFC)>logFC_filter & PValue<=p_value)

print('Filtering Result')
head(top_res_filtered)
dim(top_res_filtered)


up_regulated <- subset(top_res_filtered, logFC > 0)
down_regulated <- subset(top_res_filtered, logFC < 0)

up_genes <- up_regulated$Gene
down_genes <- down_regulated$Gene

print('up/down reg gene count')
length(up_genes)
length(down_genes)


write.table(top50_res_table, file=paste0(output_f_name, '_top50.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(top50_res_table, file=paste0(output_f_name, '_top50.csv'), quote = F, row.names = F)


write.table(top_res_table, file=paste0(output_f_name, '_unfiltered.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(top_res_table, file=paste0(output_f_name, '_unfiltered.csv'), quote = F, row.names = F)


write.table(top_res_filtered, file=paste0(output_f_name, '_filtered.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(top_res_filtered, file=paste0(output_f_name, '_filtered.csv'), quote = F, row.names = F)




write.table(up_regulated, file=paste0(output_f_name, '_B_upreg.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(up_regulated, file=paste0(output_f_name, '_B_upreg.csv'), quote = F, row.names = F)

write.table(up_genes, file=paste0(output_f_name, '_B_upreg_genes.txt'), sep='\n',quote = FALSE,
            row.names = F, col.names = F)




write.table(down_regulated, file=paste0(output_f_name, '_B_downreg.tsv'), sep='\t',quote = FALSE,
            row.names = F)
write.csv(down_regulated, file=paste0(output_f_name, '_B_downreg.csv'), quote = F, row.names = F)

write.table(down_genes, file=paste0(output_f_name, '_B_downreg_genes.txt'), sep='\n',quote = FALSE,
            row.names = F, col.names = F)



# 
# # Gene ontology analysis
# 
# 
# 
# tmp <- top_res
# 
# 
# tmp$table <- rownames_to_column(tmp$table, 'Gene')
# 
# 
# top50_gene <- rownames(top50$table)
# top50_table <- top50$table
# class(top50_table)
# 
# top50_gene
# top50_table
# top50_table_lst <- as.list(top50_table)
# top50_table_lst
# 
# # gene_name_df <- data.frame(Symbol=top50_gene)
# # gene_name_df
# # top50_table <- rownames_to_column(top50_table, 'Gene')
# # top50_table
# 
# de_lst <- list(gene_name_df, top50_table)
# de_lst
# names(de_lst) <- c('genes', 'table')
# de_lst
# 
# # rownames(top_res)
# # 
# # et_no_rep
# # 
# top50_gene
# 
# 
# 
# go_only_geneName <- goana(top50_gene)
# ?goana
# go <- goana(top50_table_lst)
# go
# 
# 
# go_only_geneName
# 
# topGO(go_only_geneName, ontology = 'BP', n=30, truncate.term = 30)
# topGO(go, ontology = 'BP', n=30, truncate.term = 30)
# 
# 

























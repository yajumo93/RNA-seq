
# root - grp_dir - files ...
# 구조의 데이터 RNAseq

getwd()

root_dir = 'C:/AMC_proj/stemcell/RNAseq/RNAseq/gdc/count/Origin_MDS'
save_path = 'table_data/mt_mut_table.tsv'

# root_dir = 'C:/AMC_proj/stemcell/RNAseq/RNAseq/gdc/count/mt_mut'
# save_path = 'table_data/mt_mut_table.tsv'

old.path <- setwd(root_dir)
getwd()


count_file_lst <- list.files(recursive = T, pattern = '*.count')
count_file_lst
count_file_lst[1]

first.sample.n = read.delim(count_file_lst[1],
                            header = F, row.names = 1)

head(first.sample.n)
class(first.sample.n) # df

count.table = data.frame(first.sample.n)
head(count.table)
str(count.table)
class(count.table)


# gsub(찾을 것, 바꿀 것, 열 지정)
count_file_lst
# name_list = gsub(".count","", count_file_lst)
# name_list

grp_vec <- c()
sam_name_vec <- c()

for (f_tag in count_file_lst){
  split_vec1 <- strsplit(f_tag, split = '/')
  org <- split_vec1[[1]][1]
  f_name <- split_vec1[[1]][2]
  f_name = gsub(".count","", f_name)
  
  print(org)
  print(f_name)
  
  grp_vec <- c(grp_vec, org)
  sam_name_vec <- c(sam_name_vec, f_name)
}

grp_vec
sam_name_vec

# library(tibble)

for (i in count_file_lst[2:length(count_file_lst)]) {
  column.n = read.delim(i, header = F, row.names = 1)
  # print(tibble(column.n))
  count.table = cbind(count.table, s = column.n)
}

colnames(count.table) = sam_name_vec

head(count.table)



# write.csv(count.table, file=file.path(save_dir, save_csv_name))
# write.table(count.table, file=file.path(save_dir, save_tsv_name), sep = "\t")
# write.table(count.table, file=file.path(save_dir, save_gsea_name), sep = "\t")
write.table(count.table, file=save_path, sep = "\t")


# DEG 분석


library(edgeR)
count.table
cpm_count <- cpm(count.table)
cpm_count
thresh <- cpm_count > 0.3 
keep <- rowSums(thresh) >= 2
keep

counts.keep <- count.table[keep,] # keep에 해당하는 행(gene)만 뽑기
head(counts.keep)

dim(count.table)
dim(counts.keep) 
head(counts.keep)
grp_vec
d0 <- DGEList(counts.keep, group = grp_vec) # colSums(counts.keep)이 lib.size임
# count_df와 sample_df로 이루어 져 있다.
head(d0)
class(d0)
str(d0)

y = calcNormFactors(d0) 
y$samples
y$samples$group

plotMDS(y)

























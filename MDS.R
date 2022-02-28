getwd()
root_dir = 'C:/AMC_proj/stemcell/RNAseq/RNAseq/gdc/count/Origin_MDS'

old.path <- setwd(root_dir)
getwd()

count_file_lst <- list.files(recursive = T)
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

for (i in count_file_lst[2:length(count_file_lst)]) {
  column.n = read.delim(i, header = F, row.names = 1)
  print(column.n)
  count.table = cbind(count.table, s = column.n)
}

colnames(count.table) = sam_name_vec

head(count.table)

save_path = 'table_data/origin_table.tsv'

# write.csv(count.table, file=file.path(save_dir, save_csv_name))
# write.table(count.table, file=file.path(save_dir, save_tsv_name), sep = "\t")
# write.table(count.table, file=file.path(save_dir, save_gsea_name), sep = "\t")
write.table(count.table, file=save_path, sep = "\t")


# plotting MDS


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

d0 <- DGEList(counts.keep, group = grp_vec) # colSums(counts.keep)이 lib.size임
# count_df와 sample_df로 이루어 져 있다.
head(d0)
class(d0)
str(d0)

y = calcNormFactors(d0) 

plotMDS(y)

setwd(old.path)



























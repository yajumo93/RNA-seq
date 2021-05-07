
## make data.frame for each htseq-count files

data_path = "D:\\junmo\\wd\\RNAseq\\data\\count_data\\RNAseq1_210419_32sample\\"
save_path = "D:\\junmo\\wd\\RNAseq\\data\\df_data\\RNAseq1_210419_32sample\\"

count_file_lst = list.files(path = data_path, 
                            pattern = "*.count")
count_file_lst[1]
first.sample.n = read.delim(paste0(data_path, count_file_lst[1]),
                            header = F, row.names = 1)

head(first.sample.n)
class(first.sample.n) # df

count.table = data.frame(first.sample.n)
head(count.table)
str(count.table)
class(count.table)


# gsub(찾을 것, 바꿀 것, 열 지정)
name_list = gsub(".count","", count_file_lst)
name_list

count_file_lst

for (i in count_file_lst[2:length(count_file_lst)]) {
  path_list = paste0(data_path, i)
  column.n = read.delim(path_list, header = F, row.names = 1)
  count.table = cbind(count.table, s = column.n)
}

colnames(count.table) = name_list

head(count.table)



# ## combine normal & tumor
# count.table.nt = cbind(count.table.n, count.table)
# #colnames(count.table.nt) = c(rep("n", 5), rep("t", 44))
# head(count.table.nt)



##save file

write.csv(count.table, paste0(save_path, "count_table.csv"))

write.table(count.table, file=paste0(save_path, "count_table.txt"), sep = "\t")
write.table(count.table, file=paste0(save_path, "count_table.gsea.txt"), sep = "\t")























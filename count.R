
## make data.frame for each htseq-count files

data_dir = "E:/stemcell/RNAseq/gdc/count/mt_mut/clone_comp/Tera"
save_dir_name = "table_data"

# save_csv_name = 'hIPS29-EB_count.csv'
# save_tsv_name = 'hIPS29-EB_count.tsv'
# save_gsea_name = 'hIPS29-EB_count.txt'
# save_csv_name = 'Tera_IPS_count.csv'
# save_tsv_name = 'Tera_IPS_count.tsv'
# save_gsea_name = 'Tera_IPS_count.txt'
save_csv_name = 'Tera_count.csv'
save_tsv_name = 'Tera_count.tsv'
save_gsea_name = 'Tera_count.txt'

# save_dir = paste0(data_dir, save_dir_name)
save_dir = file.path(data_dir, save_dir_name)
save_dir
old.path <- setwd(data_dir)
getwd()

if(!file.exists(save_dir)){
    dir.create(save_dir_name)
}

file.exists(save_dir)

count_file_lst = list.files(path = data_dir, 
                            pattern = "*.count")
count_file_lst[1]
count_file_lst
first.sample.n = read.delim(file.path(data_dir, count_file_lst[1]),
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
  path_list = list.files(data_dir, i)
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

write.csv(count.table, file=file.path(save_dir, save_csv_name))
write.table(count.table, file=file.path(save_dir, save_tsv_name), sep = "\t")
write.table(count.table, file=file.path(save_dir, save_gsea_name), sep = "\t")


setwd(old.path)


















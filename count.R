
## make data.frame for each htseq-count files

data_dir = "/data_244/RNA/mapped/count/"
save_dir_name = "table_data/"

save_dir = paste0(data_dir, save_dir_name)
save_dir
setwd(data_dir)
getwd()

if(!file.exists(save_dir)){
    dir.create(save_dir_name)
}

file.exists(save_dir)


count_file_lst = list.files(path = data_dir, 
                            pattern = "*.count")
count_file_lst[1]
count_file_lst
first.sample.n = read.delim(paste0(data_dir, count_file_lst[1]),
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
  path_list = paste0(data_dir, i)
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

write.csv(count.table, paste0(save_dir, "count_table.csv"))

write.table(count.table, file=paste0(save_dir, "count_table.txt"), sep = "\t")
write.table(count.table, file=paste0(save_dir, "count_table.gsea.txt"), sep = "\t")























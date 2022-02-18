
# GLM 방식으로 EdgeR
# https://www.youtube.com/watch?v=onGKdTg-_L4&list=PLr-cdL46Ks3U9lMIxfy3d-GnXF7qZWWcB&index=6 참조

working_dir <- "E:/stemcell_ips/gdc/RNA_seq/count/hIPS29_Teratoma/table_data"

sample_id <- c('29-A-P50-Tera', '29-B-P50-Tera')
sample_class <- c('control', 'mutation')


old.path <- setwd(working_dir)

count_data <- read.table('29-A-Teratoma_count.tsv', header = T)
count_data

head(count_data)

sam.id <- colnames(count_data)
sam.id

library(edgeR)

meta_data <- data.frame(sam.id, sample_class)

target <- meta_data

target$sam.id <- factor(target$sam.id)
target$sample_class <- factor(target$sample_class)
# levels(target$sample_class)
str(target)
target

condition <- target$sample_class
condition


# 선형 모델에 데이터 fitting
design <- model.matrix(~condition, data = target)
design

?DGEList

y <- DGEList(counts = count_data)
y <- calcNormFactors(y)
y$counts
y <- estimateDisp(y, design=design)

# setwd(old.path)























### GSEA

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("edgeR")

library(edgeR)

data_dir = "E:/stemcell/RNAseq/gdc/count/table_data/"
table_name = "Tera_IPS_count.txt"

output_f_name = 'IPS_Tera_DEG'

logFC_filter <- 1.5 # abs(logFC)
p_value <- 0.05


setwd(data_dir)
getwd()

# ?read.delim

counts <- read.delim(paste0(data_dir, table_name), row.names = 1)

raw_count <- rowSums(counts) > 0 # 최소 한쪽에서 발현이 0보다 많이 된것
raw_count
raw.over.0 <- counts[raw_count,]
dim(raw.over.0)

head(counts)
samples = colnames(counts)
sum(counts$hiPS21.B) # 87278916
sum(counts$Teratoma.13)


# grouping
class(samples)
samples[1]

grp = c()

for(s in samples){
  print(s)
  splited = strsplit(s, split=".", fixed=TRUE)
  # print(splited[[1]][1]) # 첫번째 key [[1]]의 타입은 벡터임. 첫번째 요소를 꺼낼꺼라 [1]로 꺼냄.
  
  if(splited[[1]][1] == "Teratoma"){
    grp = c(grp, 2)
  }else {
     grp = c(grp, 1)
  }
  # break
}

table(grp == 1) # origin: 1 = 21개 // Teratoma: 2 = 14개 (210629)


# class(counts) # df


#d0 <- DGEList(counts)
#y = calcNormFactors(d0)
#GSEA_table <- cpm(y, log=TRUE)

# ?cpm
# CPM : count per million
cpm_count <- cpm(counts) # 전체 fragment(전체 gene에서 발현된 모든 fragments) 백만개당 target fragment의 개수 정규화 
                         # counts에서 --> (하나의 gene / 전체 gene 합) 값임. 퍼센트는 아님
                         # 전체 발현량 : 100만 = 특정 gene 발현량 : 정규화 값(x)
head(cpm_count)

thresh <- cpm_count > 0.3 # hESO14 샘플의 경우 43개 이상이면 pass (y = (총 합 * 0.5 / 백만))
# 트리밍할 '유전자' 정해주는 과정
head(thresh)


# 테라토마나 ips 매치되는거 생각해서 동량의 것만 비교해야지. 
# 걍 이렇게만 하면 아에그냥 테라토마든 뭐든 
## 일단 cpm 역치 기준 넘는거만 전체적으로 쭉 보는거임

head(rowSums(thresh)) # 수치값: 유전자당 cpm 역치 넘는 sample 수. {gene type : 역치 넘는 sample수} 형태 (즉, max = sample 수)
# class(rowSums(thresh)) # "numeric"

table(rowSums(thresh)) # 위의 발현한 '수' 기준으로 빈도 출력. {해당되는 샘플수 : 샘플수에 해당하는 gene type 수}
# class(table(rowSums(thresh))) # table


# {발현한 수 : 그 수에 해당하는 유전자 숫자}

# rowSums(thresh) >= 2 # 유전자들에 대하여 최소 두개 샘플에서 CPM 역치 넘는 유전자만 추출
keep <- rowSums(thresh) >= 2
# keep <- rowSums(thresh) >= 1 # 유전자들에 대하여 최소 한개 샘플에서 CPM 역치 넘는 유전자만 추출

dim(counts.keep) # 숫자 일치 확인

summary(keep) # 두개 이상일 때 = {F:33890 / T:21997} // 한개 이상일 때 {F:31815 / T:24072}

head(keep)

counts.keep <- counts[keep,] # keep에 해당하는 행(gene)만 뽑기
head(counts.keep)

dim(counts)
dim(counts.keep) # cpm 역치 넘은 '유전자'만 남음.  (55887 -> 21997)

max((counts.keep)[,1]) # 첫번째 sample에서 가장 많이 발현한 값 = 14496456
max((cpm(counts.keep))[,1]) # 첫번째 sample에서 가장 많이 발현한 cpm 값 = 166241.9

# pdf("E:/stemcell/RNAseq/gdc/count/table_data/plot.pdf", width = 12, height = 12)
# plot(cpm(counts.keep)[, 1], counts.keep[,1]) # 당연히 비례 그래프 출력 
# dev.off()

# install.packages('txtplot')
# library('txtplot')
# txtplot(cpm(counts.keep)[, 1], counts.keep[,1])


#convert counts to DGEList object
# ?DGEList
# ?filterByExpr # 위 과정은 수동필터, 얘는 자동필터인듯 
# ?model.matrix
# ?estimateDisp
# ?exactTest


# group 정보를 매칭되는 샘플별로 넣어주면 되나?
# keep 자체가 cpm을 활용해 걸러낸 count table임.
d0 <- DGEList(counts.keep, group = grp) # colSums(counts.keep)이 lib.size임
# count_df와 sample_df로 이루어 져 있다.
head(d0)
class(d0)
str(d0)

head(d0$counts[, 'Teratoma.9'])
d0$samples
class(d0$samples) # df

y = calcNormFactors(d0) # normalization
head(y)
y$samples
sample_label <- c(h21.B)

plotMD(y)
plotMDS(y)


y$samples
write.table(y$samples, file='gsea.input.filtered.txt', sep='\t',quote = FALSE)
write.csv(y$samples, file = 'sample_table.csv', quote = F)

cpm_count <- cpm(y, log=T)
head(cpm_count)

write.table(cpm_count, file='gsea.input.filtered.txt', sep='\t',quote = FALSE)

# tmp_count = cpm(counts.keep, log=T) # 위에꺼랑 값이 다름
# head(tmp_count)


#########################################
# 자동필터

# head(counts)
# y = DGEList(count=counts, group=grp)
# class(y)
# y
# keep_mtd <- filterByExpr(y, group=grp)
# y = counts[keep_mtd, ]
# head(y)
# class(y)
# summry(y)
#########################################

# getwd()

y1 = estimateDisp(y) # To estimate common dispersion and tagwise dispersions in one run

y1 # count, sample에다가 추가로 분산 관련해서 더 붙음. 21992 (row)는 gene의 갯수임


?exactTest

# et = exactTest(y) # err

et = exactTest(y1)

et
plotMD(et)

top50 <- topTags(et, n = 50)
top50

plotMD(top50)


# et2 = exactTest(y1, pair = c(1, 2))

# et2

# topTags(et2)



y1$samples$group


row_len <- dim(et$table)[1]
row_len

top_res <- topTags(et, n = row_len) # fdr 값 새겨진 table get
top_res_table <- top_res$table
dim(top_res_table)
head(top_res_table)
et$comparison # 1: ips, 2: Tera

library(tibble)

top_res_table <- rownames_to_column(top_res_table, 'Gene')
# top50_res_table <- rownames_to_column(top50_table, 'Gene')

head(top_res_table)

# head(top50_res_table)
# dim(top50_res_table)



top_res_filtered <- subset(top_res_table, 
                           abs(logFC)>logFC_filter & PValue<=p_value &
                             FDR<=0.05)

dim(top_res_filtered)

getwd()

up_regulated <- subset(top_res_filtered, logFC > 0)
down_regulated <- subset(top_res_filtered, logFC < 0)

up_genes <- up_regulated$Gene
down_genes <- down_regulated$Gene

up_genes



# write.table(top50_res_table, file=paste0(output_f_name, '_top50.tsv'), sep='\t',quote = FALSE,
#             row.names = F)
# write.csv(top50_res_table, file=paste0(output_f_name, '_top50.csv'), quote = F, row.names = F)




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















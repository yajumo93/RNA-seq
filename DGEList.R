### GSEA

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# BiocManager::install("edgeR")

library(edgeR)

data_path = "D:\\junmo\\wd\\RNAseq\\data\\df_data\\RNAseq1_210419_32sample\\"
table_name = "count_table.gsea.txt"

# setwd(data_path)
# getwd()

?read.delim

counts <- read.delim(paste0(data_path, table_name), row.names = 1)
head(counts)


#d0 <- DGEList(counts)
#y = calcNormFactors(d0)
#GSEA_table <- cpm(y, log=TRUE)
?cpm
GSEA_table <- cpm(counts) # 전체 fragment 백만개당 target fragment의 개수 정규화 
head(GSEA_table)

thresh <- GSEA_table > 0.5 # 0.5%니까 5천개 이상 카운트로 역치줌
# 트리밍할 '유전자' 정해주는 과정
head(thresh)


# 테라토마나 ips 매치되는거 생각해서 동량의 것만 비교해야지. 
# 걍 이렇게만 하면 아에그냥 테라토마든 뭐든 
## 일단 cpm 역치 기준 넘는거만 전체적으로 쭉 보는거임

head(rowSums(thresh)) # 수치값: 유전자당 샘플에서 발현한 카운트 중 cpm 역치 넘는 값 (즉, max = sample 수)
table(rowSums(thresh)) # 위의 발현한 '수' 기준으로 빈도 출력. 
# {발현한 수 : 그 수에 해당하는 유전자 숫자}

rowSums(thresh) >= 2 # 샘플에서 발현한 유전자 중 최소 두개는 cpm 역치 넘는 유전자만 TRUE
# keep <- rowSums(thresh) >= 2
keep <- rowSums(thresh) >= 1

summary(keep)

counts.keep <- counts[keep,]
head(counts.keep)

dim(counts)
dim(counts.keep) # cpm 역치 넘은 '유전자'만 남음. 

max((cpm(counts.keep))[,1]) # 첫번째 sample에서 가장 많이 발현한 값

plot(cpm(counts.keep)[, 1], counts.keep[,1]) # 당연히 비례 그래프 출력 


#convert counts to DGEList object
?DGEList
?filterByExpr # 위 과정은 수동필터, 얘는 자동필터인듯 
?model.matrix
?estimateDisp
?exactTest

# group 정보를 매칭되는 샘플별로 넣어주면 되나?
d0 <- DGEList(counts.keep) # colSums(counts.keep)이 lib.size임
d0
class(d0)
str(d0)



head(d0$counts)
d0$samples

y = calcNormFactors(d0) # normalization
head(y)

GSEA_table <- cpm(y, log=T)
GSEA_table

write.table(GSEA_table, file='gsea.input.filtered.txt', sep='\t',quote = FALSE)


































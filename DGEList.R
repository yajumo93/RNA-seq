### GSEA

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# BiocManager::install("edgeR")

library(edgeR)

data_dir = "/data_244/RNA/mapped/count/table_data/"
table_name = "count_table.gsea.txt"

setwd(data_dir)
getwd()

# ?read.delim

counts <- read.delim(paste0(data_dir, table_name), row.names = 1)
head(counts)
samples = colnames(counts)
sum(counts$hESO14) # 87278916
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
head(cpm_count)

thresh <- cpm_count > 0.5 # hESO14 샘플의 경우 43개 이상이면 pass (y = (총 합 * 0.5) / 백만)
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

summary(keep) # 두개 이상일 때 = {F:33890 / T:21997} // 한개 이상일 때 {F:31815 / T:24072}

head(keep)

counts.keep <- counts[keep,] # keep에 해당하는 행(gene)만 뽑기
head(counts.keep)

dim(counts)
dim(counts.keep) # cpm 역치 넘은 '유전자'만 남음.  (55887 -> 21997)

max((counts.keep)[,1]) # 첫번째 sample에서 가장 많이 발현한 값 = 14496456
max((cpm(counts.keep))[,1]) # 첫번째 sample에서 가장 많이 발현한 cpm 값 = 166241.9

plot(cpm(counts.keep)[, 1], counts.keep[,1]) # 당연히 비례 그래프 출력 

# install.packages('txtplot')
# library('txtplot')
txtplot(cpm(counts.keep)[, 1], counts.keep[,1])


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

topTags(et)


# et2 = exactTest(y1, pair = c(1, 2))

# et2

# topTags(et2)



y1$samples$group


















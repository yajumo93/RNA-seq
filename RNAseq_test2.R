
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("edgeR")


library(edgeR)


# Reading in the data

old.path <- setwd('E:/stemcell_ips/gdc/RNA_seq/EdgeR_test_data')

rawdata <- read.delim("tabelS1_test.txt", check.names=FALSE, 
                      stringsAsFactors=FALSE)

head(rawdata)

# library(edgeR)

y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
y


# Annotation

# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
idfound
y <- y[idfound,]
dim(y)


egREFSEQ <- toTable(org.Hs.egREFSEQ)
head(egREFSEQ)
tail(egREFSEQ)

m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
head(y$genes)
str(y)



# Filtering and normalization

head(y$counts)

rowSums(y$counts)
o <- order(rowSums(y$counts), decreasing=TRUE)
o
y <- y[o,]
y
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)



y$samples$lib.size <- colSums(y$counts)
y$samples$lib.size


rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL



y <- calcNormFactors(y)
y$samples

plotMDS(y)


#  The design matrix

Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)


design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
design


# Estimating the dispersion

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion


plotBCV(y)


fit <- glmFit(y, design)
# fit <- glmQLFit(y, design, robust=T)
fit




# 4.1.8 Differential expression

lrt <- glmLRT(fit)
lrt
summary(lrt)
class(lrt$genes)
class(lrt$table)
str(lrt)
lrt$genes

topTags(lrt)
class(topTags(lrt))



colnames(design)

o <- order(lrt$table$PValue)
o

cpm(y)[o[1:10],]

summary(decideTests(lrt))


plotMD(lrt)

abline(h=c(-1, 1), col='blue')


BiocManager::install("GO.db")

go <- goana(lrt)
class(lrt)
class(lrt$genes$)

go



















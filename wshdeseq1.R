setwd("C:/Users/a/Desktop/数据分析/wsh")

library(DESeq2)
library(dplyr)
library(readxl)
#导入counts数据矩阵，以行为基因，列为样本


data2co <- read.xlsx('C:/Users/a/Desktop/数据分析/wsh/GSE185694_fpkm_2.xlsx', colNames = TRUE,rowNames = TRUE)
sdata2co <- data2co[c(grep("Con",colnames(data2co)))]
sdata2co <- sdata2co[-c(grep("FPKM",colnames(sdata2co)))]

data2fe <- read.xlsx('C:/Users/a/Desktop/数据分析/wsh/GSE185694_fpkm_2.xlsx', colNames = TRUE,rowNames = TRUE)
sdata2fe <- data2fe[c(grep("Fe",colnames(data2fe)))]
sdata2fe <- sdata2fe[-c(grep("FPKM",colnames(sdata2fe)))]

scb2data <- cbind(sdata2co, sdata2fe)

count <- scb2data
count <- round(count)
## 过滤在所有重复样本中小于1的基因，表达量太低也没研究意义
#count <- count[rowMeans(count)>1,]
##载入样本信息
fzdata2 <- read.csv('wshfe2注释.csv', sep = ',',header = T,row.names = 1)

sfzdata2 <- fzdata2[c(grep("count",fzdata2$digital)),]
#fzdata2 <- read.table("C:/Users/a/Desktop/表达矩阵/注释信息.txt",header = T,row.names = 1)

#一定要变为因子数据，否者用DESeq2包分析时候会出错
sfzdata2[,2] <- as.factor(sfzdata2$type)


all(rownames(sfzdata2) %in% colnames(count))
all(rownames(sfzdata2) == colnames(count))

dds <-  DESeqDataSetFromMatrix(countData = count,colData = sfzdata2,design = ~type)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"all_diff_fe.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diff,"all_diffsigcd4.csv")
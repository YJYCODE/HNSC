

library(data.table)
library(limma)
cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
inputFile=paste0(path,cancer,"\\1-data\\",cancer,"_FPKM.txt")
#读取表达矩阵,整理
data=as.matrix(fread(inputFile,sep="\t",header=T,check.names=F))

rownames(data)=as.character(data[,2])
exp=data[,3:ncol(data)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) ## 对同一个基因多行的数据，取均值

out=data[rowMeans(data)>0,] # 过滤表达值低的基因(所有样本中全为0的基因)

#数据转换（转录组的RNA-seq数据转换为芯片数据）
v <-voom(data, plot = F, save.plot = F) # limma包voom函数，可以把RNA-seq的数据标准化为类似芯片数据的结果
out=v$E
out=rbind(ID=colnames(out),out)
dir.create(paste0(path,cancer,"\\4-imm\\cibersort\\"))
write.table(out,file=paste0(path,cancer,"\\4-imm\\cibersort\\","uniq.symbol.txt"),
            sep="\t",quote=F,col.names=F,row.names = T)

setwd(paste0(path,cancer,"\\4-imm\\cibersort\\"))
#运行CIBERSORT，得到免疫细胞含量结果
source(paste0(path,"CIBERSORT.R"))
results=CIBERSORT(paste0(path,"ref.txt"), "uniq.symbol.txt", perm=100, QN=TRUE)
#perm置换次数=100，QN分位数归一化=TRUE
# 结果文件：
# RMSE就是均方根误差(Root Mean Squared Error)
# p值，小于0.05就是计算可信
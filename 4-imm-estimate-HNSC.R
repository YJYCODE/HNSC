


#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

#载入包
library(limma) ##avereps函数
library(estimate)
library(data.table)

cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
inputFile=paste0(path,cancer,"\\1-data\\",cancer,"_FPKM.txt")
sample_group <- read.table(paste0(path,cancer,"\\1-data\\",cancer,"_FPKM_group.txt"),
                           sep = "\t",header = T)
#读取表达矩阵,整理
data=as.matrix(fread(inputFile,sep="\t",header=T,check.names=F))

rownames(data)=as.character(data[,2])
exp=data[,3:ncol(data)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) ## 对同一个基因多行的数据，取均值

out=data[rowMeans(data)>0,] # 过滤表达值低的基因(所有样本中全为0的基因)


#输出整理后的矩阵文件
out=rbind(ID=colnames(out),out)
dir.create(paste0(path,cancer,"\\4-imm\\"))
write.table(out,file=paste0(path,cancer,"\\4-imm\\","uniq.symbol.txt"),
            sep="\t",quote=F,col.names=F)

#运行estimate包
setwd(paste0(path,cancer,"\\4-imm\\"))
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")

#输出每个样品的打分
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores <- cbind(rownames(scores),scores)
colnames(scores) <- c("sample","StromalScore","ImmuneScore","ESTIMATEScore")
out <- merge(scores,sample_group[,c(1,7)],by.x="sample",by.y="sample_barcode",
             all.x=T)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=T,row.names = F)

############## 箱线图 ###################################

##基质分数
wilcox.test(as.numeric(out[out$group=="ecDNA -",2]),
            as.numeric(out[out$group=="ecDNA +",2]))
boxplot(as.numeric(out[out$group=="ecDNA -",2]),
        as.numeric(out[out$group=="ecDNA +",2]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="StromalScore")

##免疫分数
wilcox.test(as.numeric(out[out$group=="ecDNA -",3]),
            as.numeric(out[out$group=="ecDNA +",3]))
boxplot(as.numeric(out[out$group=="ecDNA -",3]),
        as.numeric(out[out$group=="ecDNA +",3]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="ImmuneScore")

## 均无差异

##总分数
wilcox.test(as.numeric(out[out$group=="ecDNA -",4]),
            as.numeric(out[out$group=="ecDNA +",4]))
boxplot(as.numeric(out[out$group=="ecDNA -",4]),
        as.numeric(out[out$group=="ecDNA +",4]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="ESTIMATEScore")














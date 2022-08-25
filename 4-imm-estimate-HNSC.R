
library(limma) 
library(estimate)
library(data.table)

cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
inputFile=paste0(path,cancer,"\\1-data\\",cancer,"_FPKM.txt")
sample_group <- read.table(paste0(path,cancer,"\\1-data\\",cancer,"_FPKM_group.txt"),
                           sep = "\t",header = T)

data=as.matrix(fread(inputFile,sep="\t",header=T,check.names=F))

rownames(data)=as.character(data[,2])
exp=data[,3:ncol(data)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) 

out=data[rowMeans(data)>0,] 


out=rbind(ID=colnames(out),out)
dir.create(paste0(path,cancer,"\\4-imm\\"))
write.table(out,file=paste0(path,cancer,"\\4-imm\\","uniq.symbol.txt"),
            sep="\t",quote=F,col.names=F)


setwd(paste0(path,cancer,"\\4-imm\\"))
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")


scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores <- cbind(rownames(scores),scores)
colnames(scores) <- c("sample","StromalScore","ImmuneScore","ESTIMATEScore")
out <- merge(scores,sample_group[,c(1,7)],by.x="sample",by.y="sample_barcode",
             all.x=T)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=T,row.names = F)


##基质分数
wilcox.test(as.numeric(out[out$group=="ecDNA -",2]),
            as.numeric(out[out$group=="ecDNA +",2]))
boxplot(as.numeric(out[out$group=="ecDNA -",2]),
        as.numeric(out[out$group=="ecDNA +",2]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="StromalScore")


wilcox.test(as.numeric(out[out$group=="ecDNA -",3]),
            as.numeric(out[out$group=="ecDNA +",3]))
boxplot(as.numeric(out[out$group=="ecDNA -",3]),
        as.numeric(out[out$group=="ecDNA +",3]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="ImmuneScore")



wilcox.test(as.numeric(out[out$group=="ecDNA -",4]),
            as.numeric(out[out$group=="ecDNA +",4]))
boxplot(as.numeric(out[out$group=="ecDNA -",4]),
        as.numeric(out[out$group=="ecDNA +",4]),
        names=c("ecDNA -","ecDNA +"),
        col=c("blue","red"),
        pch=20,
        main="ESTIMATEScore")














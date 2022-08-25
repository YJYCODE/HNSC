
library(limma) 
library(estimate)
library(data.table)
library(vioplot)

cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
inputFile=paste0(path,cancer,"\\1-data\\",cancer,"_FPKM.txt")
sample_group <- read.table(paste0(path,cancer,"\\1-data\\",cancer,"_FPKM_group.txt"),
                           sep = "\t",header = T)

pFilter=0.05         

immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


Neg=sample_group[sample_group$group=="ecDNA -",1]
Pos=sample_group[sample_group$group=="ecDNA +",1]


NegImm <- intersect(row.names(immune),Neg)
PosImm <- intersect(row.names(immune),Pos)

rt=rbind(immune[NegImm,],immune[PosImm,])
NegNum=length(NegImm)
PosNum=length(PosImm)


setwd(paste0(path,cancer,"\\4-imm\\cibersort\\"))
outTab=data.frame()
pdf("vioplot.pdf",height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")


for(i in 1:ncol(rt)){
  if(sd(rt[1:NegNum,i])==0){
    rt[1,i]=0.001
  }
  if(sd(rt[(NegNum+1):(NegNum+PosNum),i])==0){
    rt[(NegNum+1),i]=0.001
  }
  lowData=rt[1:NegNum,i]
  highData=rt[(NegNum+1):(NegNum+PosNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()


write.table(outTab,
            file=paste0(path,cancer,"\\4-imm\\cibersort\\diff.result.txt"),
            sep="\t",row.names=F,quote=F)


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




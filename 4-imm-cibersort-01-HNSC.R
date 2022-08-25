

library(data.table)
library(limma)
cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
inputFile=paste0(path,cancer,"\\1-data\\",cancer,"_FPKM.txt")
#��ȡ�������,����
data=as.matrix(fread(inputFile,sep="\t",header=T,check.names=F))

rownames(data)=as.character(data[,2])
exp=data[,3:ncol(data)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) ## ��ͬһ��������е����ݣ�ȡ��ֵ

out=data[rowMeans(data)>0,] # ���˱���ֵ�͵Ļ���(����������ȫΪ0�Ļ���)

#����ת����ת¼���RNA-seq����ת��ΪоƬ���ݣ�
v <-voom(data, plot = F, save.plot = F) # limma��voom���������԰�RNA-seq�����ݱ�׼��Ϊ����оƬ���ݵĽ��
out=v$E
out=rbind(ID=colnames(out),out)
dir.create(paste0(path,cancer,"\\4-imm\\cibersort\\"))
write.table(out,file=paste0(path,cancer,"\\4-imm\\cibersort\\","uniq.symbol.txt"),
            sep="\t",quote=F,col.names=F,row.names = T)

setwd(paste0(path,cancer,"\\4-imm\\cibersort\\"))
#����CIBERSORT���õ�����ϸ���������
source(paste0(path,"CIBERSORT.R"))
results=CIBERSORT(paste0(path,"ref.txt"), "uniq.symbol.txt", perm=100, QN=TRUE)
#perm�û�����=100��QN��λ����һ��=TRUE
# ����ļ���
# RMSE���Ǿ��������(Root Mean Squared Error)
# pֵ��С��0.05���Ǽ������
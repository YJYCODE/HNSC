
library(data.table)
path <- "E:\\Training\\12-topic-data-analysis\\"
setwd(path)
cancer <- "HNSC"
FPKM <- as.data.frame(fread(paste0(path,cancer,"\\","1-data\\",
                                   cancer,"_FPKM.txt"),header = T,sep = "\t"))


FPKM <- FPKM[rowMeans(FPKM[,3:ncol(FPKM)],na.rm = T)>0,] 
FPKM$mean_negative <- apply(FPKM[,3:115],1,function(x){mean(2^x-1,na.rm = T)})
FPKM$mean_positive <- apply(FPKM[,116:154],1,function(x){mean(2^x-1,na.rm = T)})
FPKM$logFC <- log2(FPKM$mean_positive)-log2(FPKM$mean_negative)
FPKM$P_value <- apply(FPKM[,3:154],1,
                      function(x) wilcox.test(x[1:113],x[114:152])$p.value)
FPKM$FDR <- p.adjust(FPKM$P_value, method = "fdr")

out <- FPKM[,c(1,2,155:ncol(FPKM))]
out_diff <- out[out$P_value<0.05,]
out_diff <- out_diff[complete.cases(out_diff[,6]),]

dir.create(paste0(path,cancer,"\\2-expression_diff\\"))
setwd(paste0(path,cancer,"\\2-expression_diff\\"))
write.table(out,file="exp_all.txt",sep="\t",row.names=F,quote=F)
write.table(out_diff,file="exp_diff.txt",sep="\t",row.names=F,quote=F)

diff_id <- out_diff[,c(1,2,5)]
dir.create(paste0(path,cancer,"\\3-exp_GO_KEGG\\"))
setwd(paste0(path,cancer,"\\3-exp_GO_KEGG\\"))
write.table(diff_id,file="id_diff.txt",sep="\t",row.names=F,quote=F)



library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
cancer <- "HNSC"
pvalueFilter=0.05         
qvalueFilter=0.05         
path <- "E:\\Training\\12-topic-data-analysis\\"

rt=read.table(paste0(path,cancer,"\\3-exp_GO_KEGG\\","id_diff.txt"),
              sep="\t",header=T,check.names=F)
library(stringr)
splitEnsembl <- function(Ensembl){
  return(str_split(Ensembl,'[.]',simplify = T)[1]) 
}
rt$Ensembl <- sapply(rt$id,splitEnsembl,simplify = T)

library(AnnotationDbi)
library(org.Hs.eg.db)
rt$EntrezID1 <- mapIds(org.Hs.eg.db,keys = rt$name,
                      column = 'ENTREZID',
                      keytype = 'SYMBOL',
                      multiVals = 'first')
rt$EntrezID2 <- mapIds(org.Hs.eg.db,keys = rt$Ensembl,
                      column = 'ENTREZID',
                      keytype = 'ENSEMBL',
                      multiVals = 'first')
rt$EntrezID <- ifelse(is.na(rt$EntrezID1),rt$EntrezID2,rt$EntrezID1)



rt=rt[is.na(rt[,"EntrezID"])==F,]                            
gene=rt$EntrezID
geneFC=2^rt$logFC
names(geneFC)=gene

setwd(paste0(path,cancer,"\\3-exp_GO_KEGG\\"))

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, 
            pvalueCutoff =1, qvalueCutoff = 1, 
            ont="all", readable =T) 
GO=as.data.frame(kk)
GO=GO[GO$pvalue<pvalueFilter,]

write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F,col.names = T)


pdf(file="barplot.pdf",width = 9,height = 7)
bar=barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY",
            color = "pvalue") + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

GO_BP <- GO[1:10,]
library(ggplot2)
ggplot(data=GO_BP, aes(x=Count, y=Description,fill=pvalue)) +
  geom_bar(stat="identity", width=0.8) + 
  theme_test() + 
  xlab("gene number") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "Top 10 GO Terms")


library(dplyr)
library(forcats)
GO_BP %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(x=Description, y=Count,fill=pvalue)) +
  geom_bar(stat="identity",alpha=.6, width=.4) +
  theme_test() + 
  xlab("gene number") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "KEGG pathway")+
  coord_flip()


pdf(file="bubble.pdf",width = 9,height = 7)
bub=dotplot(kk,showCategory = 10, orderBy = "GeneRatio",split="ONTOLOGY", 
            color = "pvalue") + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

ggplot(GO_BP,aes(pvalue,Description))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_color_gradient(low="green",high = "red")+ #??????????
  labs(color=expression(-log[10](pvalue)),size="Count",  
       x="pvalue",y="GO name",title="Top 10 GO terms")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))




kk2 <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk2)
#KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[KEGG$pvalue<pvalueFilter,]

write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)


pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk2, drop = TRUE, showCategory = 10, color = "pvalue")
dev.off()

library(ggplot2)
ggplot(data=KEGG, aes(x=Count, y=Description,fill=pvalue)) +
  geom_bar(stat="identity", width=0.8) + 
  theme_test() + 
  xlab("gene number") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "KEGG pathway")


library(dplyr)
library(forcats)
KEGG %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(x=Description, y=Count,fill=pvalue)) +
  geom_bar(stat="identity",alpha=.6, width=.4) +
  theme_test() + 
  xlab("gene number") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "KEGG pathway")+
  coord_flip()


pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 10, orderBy = "GeneRatio",color = "pvalue")
dev.off()

ggplot(KEGG,aes(pvalue,Description))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_color_gradient(low="green",high = "red")+ #??????????
  labs(color=expression(-log[10](pvalue)),size="Count",  
       x="pvalue",y="Pathway name",title="KEGG Pathway enrichment")+
  theme_bw()
  

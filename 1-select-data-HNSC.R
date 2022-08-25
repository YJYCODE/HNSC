

### data select and merge ######################################################
library(data.table)
library(openxlsx)
path <- "E:\\Training\\12-topic-data-analysis\\"
setwd(path)

cancer <- "HNSC"
tag <- "Head and Neck"
FPKM_file <- paste0(path,cancer,"\\","TCGA-",cancer,".htseq_fpkm.tsv")
SNV_file <- paste0(path,cancer,"\\","TCGA-",cancer,".mutect2_snv.tsv")
survival_file <- paste0(path,cancer,"\\",cancer,"_survival.txt")

## load data
ecDNA <- read.xlsx("TCGA_ecDNA_classification.xlsx",sheet = 1)
ecDNA_cancer <- ecDNA[ecDNA$lineage==tag,]
FPKM_data <- as.data.frame(data.table::fread(FPKM_file,header=T))
SNV_data <- as.data.frame(data.table::fread(SNV_file,header=T))
survival_data <- as.data.frame(data.table::fread(survival_file,header=T))

## change TCGA code ("01A" to "01")
colnames(FPKM_data) <- substr(colnames(FPKM_data),1,15)
SNV_data$Sample_ID <- substr(SNV_data$Sample_ID,1,15)

## merge data (only tumor)
ecDNA_cancer_tumor <- ecDNA_cancer[ecDNA_cancer$tumor_or_normal=="tumor",]
test1 <- intersect(ecDNA_cancer_tumor$sample_barcode,colnames(FPKM_data))
test2 <- intersect(ecDNA_cancer_tumor$sample_barcode,SNV_data$Sample_ID)
test3 <- intersect(ecDNA_cancer_tumor$sample_barcode,survival_data$sample)

FPKM_data_use <- FPKM_data[,c("Ensembl_ID",test1)]
SNV_data_use <- SNV_data[SNV_data$Sample_ID%in%test2,]
survival_data_use <- survival_data[survival_data$sample%in%test3,]


FPKM_group <- ecDNA_cancer_tumor[ecDNA_cancer_tumor$sample_barcode%in%test1,]
FPKM_group$group <- ifelse(FPKM_group$sample_classification=="Circular",
                           "ecDNA +","ecDNA -")
table(FPKM_group$group)
# SNV和survival的分组跟FPKM一样
SNV_group <- ecDNA_cancer_tumor[ecDNA_cancer_tumor$sample_barcode%in%test2,]
SNV_group$group <- ifelse(SNV_group$sample_classification=="Circular",
                          "ecDNA +","ecDNA -")
table(SNV_group$group)
survial_group <- ecDNA_cancer_tumor[ecDNA_cancer_tumor$sample_barcode%in%test3,]
survial_group$group <- ifelse(survial_group$sample_classification=="Circular",
                              "ecDNA +","ecDNA -")
table(survial_group$group)

#### 输出数据 ###################################################
# 1. 表达矩阵，整理成:
#列：symbol,id,sample1,sample2.....，ecDNA-样本在前，ecDNA+样本在后
#行：基因
sample <- c(FPKM_group[FPKM_group$group=="ecDNA -",1],
            FPKM_group[FPKM_group$group=="ecDNA +",1])
FPKM_data_final <- FPKM_data_use[,c("Ensembl_ID",sample)]
# 基因id注释
gtf <- as.data.frame(fread("gencode.v22.annotation.gtf",sep = "\t",header = F))
gene_gtf <- gtf[gtf$V3=="gene",]
rownames(gene_gtf) <- NULL
split_infor_id <- function(data){
  i <- gsub("gene_id \"(.*?)\"; gene_type \"(.*?)\"; gene_status \"(.*?)\"; gene_name \"(.*?)\";.*","\\1",data)
  return(i)
}
split_infor_name <- function(data){
  i <- gsub("gene_id \"(.*?)\"; gene_type \"(.*?)\"; gene_status \"(.*?)\"; gene_name \"(.*?)\";.*","\\4",data)
  return(i)
}
gene_gtf$id <- sapply(gene_gtf$V9,split_infor_id)
gene_gtf$name <- sapply(gene_gtf$V9,split_infor_name)

need_gtf <- gene_gtf[,10:11] 

FPKM_data_final <- merge(need_gtf,FPKM_data_final,
                         by.x = "id",by.y = "Ensembl_ID",all = T)
## 1-2列，gene name  3:115列 ecDNA- 116:154列 ecDNA+

# 2.survival数据
survival_data_final <- merge(survival_data_use,survial_group,
                             by.x = "sample",by.y = "sample_barcode",all = T)

dir.create(paste0(path,cancer,"\\","1-data\\"))
setwd(paste0(path,cancer,"\\","1-data\\"))
write.table(FPKM_data_final,file = paste0(cancer,"_FPKM.txt"),
            quote = F,row.names = F,sep = "\t",col.names = T)
write.table(FPKM_group,file = paste0(cancer,"_FPKM_group.txt"),
            quote = F,row.names = F,sep = "\t",col.names = T)
write.table(SNV_data_use,file = paste0(cancer,"_SNV.txt"),
            quote = F,row.names = F,sep = "\t",col.names = T)
write.table(SNV_group,file = paste0(cancer,"_SNV_group.txt"),
            quote = F,row.names = F,sep = "\t",col.names = T)
write.table(survival_data_final,file = paste0(cancer,"_survival.txt"),
            quote = F,row.names = F,sep = "\t",col.names = T)









































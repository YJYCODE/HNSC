

library(survival)
library(survminer)
cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
sur_file <- paste0(path,cancer,"\\1-data\\",cancer,"_survival.txt")
sur <- read.table(sur_file,header = T,sep = "\t")
sur$group <- factor(sur$group,levels = c("ecDNA -","ecDNA +"))
## 去除掉OS time 为NA的
sur_reNA <- sur[complete.cases(sur[,4]),]
fit <- survfit(Surv(OS.time, OS) ~ group, data = sur_reNA)
dir.create(paste0(path,cancer,"\\5-sur\\"))
setwd(paste0(path,cancer,"\\5-sur\\"))
png(paste0(cancer,"_OS.png"),height = 600,width = 600)
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = F,# 置信区间
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", #  添加中位生存时间线
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("blue", "red"),
           xlab = "Follow up time(days)"# 指定x轴标签
)
dev.off()

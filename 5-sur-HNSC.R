

library(survival)
library(survminer)
cancer <- "HNSC"
path <- "E:\\Training\\12-topic-data-analysis\\"
sur_file <- paste0(path,cancer,"\\1-data\\",cancer,"_survival.txt")
sur <- read.table(sur_file,header = T,sep = "\t")
sur$group <- factor(sur$group,levels = c("ecDNA -","ecDNA +"))

sur_reNA <- sur[complete.cases(sur[,4]),]
fit <- survfit(Surv(OS.time, OS) ~ group, data = sur_reNA)
dir.create(paste0(path,cancer,"\\5-sur\\"))
setwd(paste0(path,cancer,"\\5-sur\\"))
png(paste0(cancer,"_OS.png"),height = 600,width = 600)
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = F,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(), 
           palette = c("blue", "red"),
           xlab = "Follow up time(days)"
)
dev.off()

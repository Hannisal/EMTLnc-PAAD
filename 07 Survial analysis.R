# time:2021/4/7
rm(list = ls())
options(stringsAsFactors = F)
library(survival) 
library(survminer) 
inputFile="TCGA_lasso_cox.csv";riskKMFile="tcga.riskKM.pdf"#TCGA
inputFile="ICGC_lasso_cox.csv";riskKMFile="icgc.riskKM.pdf"#ICGC
rt<-read.csv(inputFile,header = T,row.names = 1,check.names = F)
## 生存曲线
fit <- survfit(Surv(OStime, OSstatus)~risk, rt)
ggsurvplot(fit,rt,risk.table=TRUE,
           conf.int=TRUE,
           palette = c("#FF1E24","#0099FF"),
           pval=TRUE,
           pval.method=TRUE)
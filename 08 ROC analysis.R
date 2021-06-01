# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
library(timeROC)
library(survival)
library(tidyverse)
riskFile="TCGA_lasso_cox.csv"
riskFile="ICGC_lasso_cox.csv"
exp<-read.csv(riskFile,header = T,row.names = 1,check.names = F)
exp_ROC<-timeROC(exp$OStime/365,exp$OSstatus,exp$riskScore,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
plot(exp_ROC,time = 1,title = F,col="#FF2216",lwd=3)
plot(exp_ROC,time = 2,title = F,col="#4F8A6C",lwd=3,add = T)
plot(exp_ROC,time = 3,title = F,col="#006499",lwd=3,add = T)
legend("bottomright",
       c(paste0("1 year AUC = ",signif(exp_ROC$AUC[1],2)),
         paste0("2 year AUC = ",signif(exp_ROC$AUC[2],2)),
         paste0("3 year AUC = ",signif(exp_ROC$AUC[3],2))),
       col = c("#FF2216","#4F8A6C","#006499"),
       bty = "n",
       lwd=3,
       cex = 1)

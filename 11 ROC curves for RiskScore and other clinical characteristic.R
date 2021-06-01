# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
# TCGA
library(timeROC)
library(survival)
library(tidyverse)
load("TCGA-riskscore-clincial.Rdata")
exp_Age<-timeROC(rt$OStime/365,rt$OSstatus,rt$Age,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_Gender<-timeROC(rt$OStime/365,rt$OSstatus,rt$Gender,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)

exp_Grade<-timeROC(rt$OStime/365,rt$OSstatus,rt$Grade,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)

exp_T<-timeROC(rt$OStime/365,rt$OSstatus,rt$T_stage,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_N<-timeROC(rt$OStime/365,rt$OSstatus,rt$N_stage,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_ROC<-timeROC(rt$OStime/365,rt$OSstatus,rt$`Risk score`,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
i=1
for (i in 1:3) {
  pdf(paste0("Multivariate ROC of",i,"Year.pdf"))
  plot(exp_ROC,time = i,title = F,col="#FF2216",lwd=3)
  plot(exp_Age,time = i,title = F,col="#267428",lwd=3,add = T)
  plot(exp_Gender,time = i,title = F,col="#2476E8",lwd=3,add = T)
  #exp_T
  plot(exp_T,time = i,title = F,col="#FFA41E",lwd=3,add = T)
  
  #exp_N
  plot(exp_N,time = i,title = F,col="#96476F",lwd=3,add = T)
  plot(exp_Grade,time = i,title = F,col="#FF92C8",lwd=3,add = T)
  legend("bottomright",
         c(paste0("Riskscore AUC = ",signif(exp_ROC$AUC[i],2)),
           paste0("Age AUC = ",signif(exp_Age$AUC[i],2)),
           paste0("Gender AUC = ",signif(exp_Gender$AUC[i],2)),
           paste0("Grade AUC = ",signif(exp_Grade$AUC[i],2)),
           paste0("T AUC = ",signif(exp_T$AUC[i],2)),
           paste0("N AUC = ",signif(exp_N$AUC[i],2))),
         col = c("#FF2216","#267428","#2476E8","#FF92C8","#FFA41E","#96476F"),
         bty = "n",
         lwd=3,
         cex = 1.2)
  dev.off()
}
# ICGC
rm(list = ls())
library(timeROC)
library(survival)
library(tidyverse)
load("ICGC-riskscore-clincial.Rdata")
exp_Age<-timeROC(rt$OStime/365,rt$OSstatus,rt$Age,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_Gender<-timeROC(rt$OStime/365,rt$OSstatus,rt$Gender,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_T<-timeROC(rt$OStime/365,rt$OSstatus,rt$T_stage,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_N<-timeROC(rt$OStime/365,rt$OSstatus,rt$N_stage,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
exp_ROC<-timeROC(rt$OStime/365,rt$OSstatus,rt$`Risk score`,cause = 1,weighting = "marginal",times = c(1,2,3),ROC = T)
for (i in 1:3) {
  pdf(paste0("ICGC Multivariate ROC of ",i," Year.pdf"))
  plot(exp_ROC,time = i,title = F,col="#FF2216",lwd=3)
  plot(exp_Age,time = i,title = F,col="#267428",lwd=3,add = T)
  plot(exp_Gender,time = i,title = F,col="#2476E8",lwd=3,add = T)
  #exp_T
  plot(exp_T,time = i,title = F,col="#FFA41E",lwd=3,add = T)
  
  #exp_N
  plot(exp_N,time = i,title = F,col="#96476F",lwd=3,add = T)
  legend("bottomright",
         c(paste0("Riskscore AUC = ",signif(exp_ROC$AUC[i],2)),
           paste0("Age AUC = ",signif(exp_Age$AUC[i],2)),
           paste0("Gender AUC = ",signif(exp_Gender$AUC[i],2)),
           paste0("T AUC = ",signif(exp_T$AUC[i],2)),
           paste0("N AUC = ",signif(exp_N$AUC[i],2))),
         col = c("#FF2216","#267428","#2476E8","#FF92C8","#FFA41E"),
         bty = "n",
         lwd=3,
         cex = 1.2)
  dev.off()
}

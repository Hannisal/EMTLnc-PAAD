# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
library(survival) 
library(survminer) 
library(pheatmap)
# The distribution of risk score, survival state, and expression heatmap of the selected lncRNAs in the training set (A) and testing set (B)
inputFile="TCGA_lasso_cox.csv";riskScoreFile="tcga.riskScore.pdf";survStatFile="tcga.survStat.pdf"
inputFile="ICGC_lasso_cox.csv";riskScoreFile="icgc.riskScore.pdf";survStatFile="icgc.survStat.pdf"
# risk score
rt<-read.csv(inputFile,header = T,row.names = 1,check.names = F)
fit <- survfit(Surv(OStime, OSstatus)~risk, rt)
rt=rt[order(rt$riskScore),]
riskClass=rt$risk
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file=riskScoreFile,width = 7,height = 3.5)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("#0099FF",lowLength),rep("#FF1E24",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("#FF1E24","#0099FF"),cex=1.2)
dev.off()
# survival state
color=as.vector(rt$OSstatus)
color[color==1]="#FF1E24"
color[color==0]="#0099FF"
pdf(file=survStatFile,width = 7,height = 3.5)
plot(rt$OStime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (Days)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#FF1E24","#0099FF"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

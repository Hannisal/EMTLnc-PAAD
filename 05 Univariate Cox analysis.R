# time:2021/4/7
rm(list = ls())
options(stringsAsFactors = F)

# TCGA
load("TCGA_exp_lnc.Rdata")
load("TCGA_clinical_Match.Rdata")

EMTlncP<-read.csv("lightcyan.csv",check.names = F,row.names = 1)
EMTlncN<-read.csv("turquoise.csv",check.names = F,row.names = 1)
EMTlnc<-c(rownames(EMTlncP),rownames(EMTlncN))

exp_EMTlnc<-exp_lnc[rownames(exp_lnc)%in%EMTlnc,]

OS<-clinical
OS<-clinical[,c(1,11,14)]
colnames(OS)<-c("ID","OSstatus","OStime")
OS$OSstatus<-ifelse(OS$OSstatus=="Dead",1,0)
OS<-OS[match(colnames(exp_EMTlnc),OS$ID),]
exp_EMTlnc<-t(exp_EMTlnc)
unicox<-cbind(OS,exp_EMTlnc)

# unicox
univar_out <- data.frame(matrix(NA,ncol(unicox)-3,6))
univar_out[,1]<-colnames(unicox)[-(1:3)]
colnames(univar_out)<-c("Genesymbol","Coeffcient","HR","lower .95","upper .95","P-value")


library(survival)
for (i in 1:(ncol(unicox)-3)){
  cox = coxph(Surv(OStime,OSstatus) ~ unicox[,i+3],data = unicox)
  cox_summ = summary(cox)
  univar_out[i,2]=cox_summ$coefficients[,1]
  univar_out[i,3]=cox_summ$coefficients[,2]
  univar_out[i,4]=cox_summ$conf.int[,3]
  univar_out[i,5]=cox_summ$conf.int[,4]
  univar_out[i,6]=cox_summ$coefficients[,5]
}
save(univar_out,file = "TCGA_univar_out.RData")
univar_out_0.05 = univar_out[univar_out[,6]<0.05,]
write.table(univar_out_0.05,file = "TCGA_univar_out_0.05.txt",sep = "\t",quote = F)
write.csv(univar_out_0.05,file = "TCGA_univar_out_0.05.csv")


# ICGC
rm(list = ls())
load("ICGC_exp_lnc.Rdata")
load("clinical.Rdata")

OS<-clinical
colnames(OS)<-c("ID","OSstatus","OStime")
OS<-OS[!OS$OStime<30,]

EMTlncP<-read.csv("lightcyan.csv",check.names = F,row.names = 1)
EMTlncN<-read.csv("turquoise.csv",check.names = F,row.names = 1)
EMTlnc<-c(rownames(EMTlncP),rownames(EMTlncN))

exp_EMTlnc<-exp_lnc[rownames(exp_lnc)%in%EMTlnc,]

exp_EMTlnc<-t(exp_EMTlnc)
rownames(exp_EMTlnc)<-pheno$icgc_donor_id
exp_EMTlnc<-exp_EMTlnc[OS$ID,]
unicox<-cbind(OS,exp_EMTlnc)



# unicox
univar_out <- data.frame(matrix(NA,ncol(unicox)-3,6))
univar_out[,1]<-colnames(unicox)[-(1:3)]
colnames(univar_out)<-c("Genesymbol","Coeffcient","HR","lower .95","upper .95","P-value")


library(survival)
for (i in 1:(ncol(unicox)-3)){
  cox = coxph(Surv(OStime,OSstatus) ~ unicox[,i+3],data = unicox)
  cox_summ = summary(cox)
  univar_out[i,2]=cox_summ$coefficients[,1]
  univar_out[i,3]=cox_summ$coefficients[,2]
  univar_out[i,4]=cox_summ$conf.int[,3]
  univar_out[i,5]=cox_summ$conf.int[,4]
  univar_out[i,6]=cox_summ$coefficients[,5]
}
save(univar_out,file = "ICGC_univar_out.RData")
univar_out_0.05 = univar_out[univar_out[,6]<0.05,]
write.table(univar_out_0.05,file = "ICGC_univar_out_0.05.txt",sep = "\t",quote = F)
write.csv(univar_out_0.05,file = "ICGC_univar_out_0.05.csv")

# venn
TCGA_unicox<-read.csv("TCGA_univar_out_0.05.csv",header = T,row.names = 1,check.names = F)
ICGC_unicox<-read.csv("ICGC_univar_out_0.05.csv",header = T,row.names = 1,check.names = F)

TCGA_unicox_risk<-TCGA_unicox[TCGA_unicox$HR>1,]
TCGA_unicox_protect<-TCGA_unicox[TCGA_unicox$HR<1,]

ICGC_unicox_risk<-ICGC_unicox[ICGC_unicox$HR>1,]
ICGC_unicox_protect<-ICGC_unicox[ICGC_unicox$HR<1,]

library(VennDiagram)
library(RColorBrewer)
venn.diagram(list(A=TCGA_unicox_risk$Genesymbol,B=ICGC_unicox_risk$Genesymbol),fill=c(brewer.pal(7,"Set1")[1:2]),
             alpha=c(0.5,0.5),cex=3,
             cat.cex=3,cat.fontface=4,lty=2,
             resolution = 300,filename = "Risklnc.tiff")
venn.diagram(list(A=TCGA_unicox_protect$Genesymbol,B=ICGC_unicox_protect$Genesymbol),fill=c(brewer.pal(7,"Set1")[1:2]),
             alpha=c(0.5,0.5),cex=3,
             cat.cex=3,cat.fontface=4,lty=2,
             resolution = 300,filename = "Protectlnc.tiff")

# Risklnc
Risklnc<-TCGA_unicox_risk$Genesymbol[which(TCGA_unicox_risk$Genesymbol%in%ICGC_unicox_risk$Genesymbol)]
# Protectlnc
Protectlnc<-TCGA_unicox_protect$Genesymbol[which(TCGA_unicox_protect$Genesymbol%in%ICGC_unicox_protect$Genesymbol)]

TCGA_unicox_venn<-TCGA_unicox[TCGA_unicox$Genesymbol%in%c(Risklnc,Protectlnc),]
ICGC_unicox_venn<-ICGC_unicox[ICGC_unicox$Genesymbol%in%c(Risklnc,Protectlnc),]

write.csv(TCGA_unicox_venn,file = "TCGA_unicox_venn.csv",row.names = F)
write.csv(ICGC_unicox_venn,file = "ICGC_unicox_venn.csv",row.names = F)

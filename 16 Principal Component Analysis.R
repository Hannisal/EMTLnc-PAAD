# time:2021/5/1
rm(list = ls())
options(stringsAsFactors = F)
# TCGA
library(tinyplanet)
library(ggplot2)
exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1)
load("TCGA_exp_mRNA.Rdata")
load("hallmark.gs.RData")

exp_emt<-exp_mRNA[rownames(exp_mRNA)%in%h$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
exp_emt<-as.matrix(exp_emt[,rownames(exp)])
group_list<-as.factor(exp$risk)
levels(group_list)<-c("low","high")

## 200 EMT related genes unicox
load("TCGA_clinical_Match.Rdata")
OS<-clinical
OS<-clinical[,c(1,11,14)]
colnames(OS)<-c("ID","OSstatus","OStime")
OS$OSstatus<-ifelse(OS$OSstatus=="Dead",1,0)
OS<-OS[match(colnames(exp_emt),OS$ID),]
exp_emt<-t(exp_emt)
unicox<-cbind(OS,exp_emt)

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
save(univar_out,file = "TCGA_EMT_univar_out.RData")
univar_out_0.05 = univar_out[univar_out[,6]<0.05,]
write.csv(univar_out_0.05,file = "TCGA_EMT_univar_out_0.05.csv")
## ICGC
rm(list = ls())
exp<-read.csv("ICGC_lasso_cox.csv",header = T,row.names = 1)
load("ICGC_exp_mRNA.Rdata")
load("hallmark.gs.RData")
load("clinical.Rdata")

colnames(exp_mRNA)<-pheno$icgc_donor_id[match(colnames(exp_mRNA),pheno$submitted_donor_id)]
exp_emt<-exp_mRNA[rownames(exp_mRNA)%in%h$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
exp_emt<-as.matrix(exp_emt[,rownames(exp)])
group_list<-as.factor(exp$risk)
levels(group_list)<-c("low","high")

## 200 EMT related genes unicox
OS<-clinical
colnames(OS)<-c("ID","OSstatus","OStime")
OS<-OS[!OS$OStime<30,]

OS<-OS[match(colnames(exp_emt),OS$ID),]
exp_emt<-t(exp_emt)
unicox<-cbind(OS,exp_emt)

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
save(univar_out,file = "ICGC_EMT_univar_out.RData")
univar_out_0.05 = univar_out[univar_out[,6]<0.05,]
write.csv(univar_out_0.05,file = "ICGC_EMT_univar_out_0.05.csv")

## prognostic EMT related genes in TCGA and ICGC
load("TCGA_EMT_univar_out.RData")
TCGA_univar_out_0.05<-univar_out[univar_out[,6]<0.05,]
load("ICGC_EMT_univar_out.RData")
ICGC_univar_out_0.05<-univar_out[univar_out[,6]<0.05,]
EMTgene<-TCGA_univar_out_0.05$Genesymbol[TCGA_univar_out_0.05$Genesymbol%in%ICGC_univar_out_0.05$Genesymbol]
## 26gene
# TCGA PCA
exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1)
load("TCGA_exp_mRNA.Rdata")
load("hallmark.gs.RData")

exp_emt<-exp_mRNA[rownames(exp_mRNA)%in%h$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
exp_emt<-as.matrix(exp_emt[,rownames(exp)])
exp_emt<-exp_emt[EMTgene,]
group_list<-as.factor(exp$risk)
levels(group_list)<-c("low","high")




pca<-prcomp(t(exp_emt),scale=T)
summary(pca)

screeplot(pca,type = "line")

# 3D
data3_center_scale <- scale(t(exp_emt), center=T, scale=T)
data3_center_scale_cov <- cov(data3_center_scale)
data3_center_scale_cov_eigen <- eigen(data3_center_scale_cov)
data3_center_scale_cov_eigen$values
data3_center_scale_cov_eigen$vectors
pc_select = 3
label = paste0("PC",c(1:pc_select))
data3_new <- data3_center_scale %*% data3_center_scale_cov_eigen$vectors[,1:pc_select]
colnames(data3_new) <- label

library(scatterplot3d)
colorl <- c("#0099FF", "#FF1E24")
# Extract same number of colors as the Group and same Group would have same color.
colors <- colorl[as.numeric(group_list)]

# 1 row 2 columns
par(mfrow=c(1,1))

scatterplot3d(data3_new, color=colors, angle=130, pch=16, main="TCGA")
legend("top", legend=levels(group_list), col=colorl, pch=16, xpd=T, horiz=T)

pca.plot = draw_pca(exp_emt,group_list)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("PCA1")+ylab("PCA2")+
  ggtitle("TCGA");pca.plot

# ICGC
rm(list = ls())
library(tinyplanet)
exp<-read.csv("ICGC_lasso_cox.csv",header = T,row.names = 1)
load("ICGC_exp_mRNA.Rdata")
load("clinical.Rdata")
exp_mRNA<-exp_mRNA[,colnames(exp_mRNA)%in%pheno$submitted_donor_id]
exp_mRNA<-exp_mRNA[,pheno$submitted_donor_id]
colnames(exp_mRNA)<-pheno$icgc_donor_id

exp_mRNA<-exp_mRNA[,rownames(exp)]
load("hallmark.gs.RData")

exp_emt<-exp_mRNA[rownames(exp_mRNA)%in%h$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
exp_emt<-as.matrix(exp_emt[,rownames(exp)])
group_list<-as.factor(exp$risk)
levels(group_list)<-c("low","high")
pca.plot = draw_pca(exp_emt,group_list)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("PCA1")+ylab("PCA2")+
  ggtitle("ICGC");pca.plot


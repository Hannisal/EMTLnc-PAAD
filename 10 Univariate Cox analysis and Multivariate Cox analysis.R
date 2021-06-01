# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
# TCGA
load("TCGA_clinical_Match.Rdata")
cl<-clinical[,c(1,6,8,9,10)]

head(cl)

exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1,check.names = F)
cl<-cl[cl$bcr_patient_barcode%in%rownames(exp),]
sum(duplicated(cl$`Patient ID`))
cl<-cl[!duplicated(cl$bcr_patient_barcode),]
cl<-cl[match(rownames(exp),cl$bcr_patient_barcode),]

cl<-cl[str_detect(cl$stage_event,"T")&str_detect(cl$stage_event,"N")&str_detect(cl$stage_event,"M"),]
cl<-cbind(cl,do.call(rbind,stringr::str_split(cl$stage_event,pattern = "T|N|M"))[,-1])
exp<-exp[rownames(exp)%in%cl$bcr_patient_barcode,]
cl<-cl[match(rownames(exp),cl$bcr_patient_barcode),]

rt<-cbind(exp[,c(1,2,(ncol(exp)-1))],cl[,-c(1,5)])
colnames(rt)[3:ncol(rt)]<-c("Risk score","Gender","Age","Grade","T_stage","N_stage","M_stage")

table(rt$Grade)
rt<-rt[!rt$Grade=="GX",]
rt$Grade<-substr(rt$Grade,2,2)
rt$Grade<-as.numeric(rt$Grade)

table(rt$T_stage)
rt<-rt[!rt$T_stage=="X",]
rt$T_stage<-as.numeric(rt$T_stage)
table(rt$N_stage)
rt<-rt[!rt$N_stage=="X",]
rt$N_stage<-substr(rt$N_stage,1,1)
rt$N_stage<-as.numeric(rt$N_stage)

table(rt$M_stage)
rt<-rt[,-ncol(rt)]

table(rt$Gender)
rt$Gender<-ifelse(rt$Gender=="FEMALE",0,1)
rt$Gender<-as.numeric(rt$Gender)
rt$Age<-as.numeric(rt$Age)
uniTab=data.frame()
library(survival)
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(OStime, OSstatus) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)

multiCox=coxph(Surv(OStime, OSstatus) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
rownames(multiTab)[1]<-"Risk score"
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="tcga.multiCox.txt",sep="\t",row.names=F,quote=F)
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  
  pdf(file=forestFile, width = 6.3,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#32804C",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
bioForest(coxFile="tcga.uniCox.txt",forestFile="tcga.uniForest.pdf",forestCol="#0D456B")
bioForest(coxFile="tcga.multiCox.txt",forestFile="tcga.multiForest.pdf",forestCol="#CC2216")
save(rt,file = "TCGA-riskscore-clincial.Rdata")
# ICGC
rm(list = ls())
cl<-read.table("donor.PACA-AU.tsv",header = T,sep = "\t")# ICGC
cl<-cl[,c(1,5,6,7,9)]
load("clinical.Rdata")
pheno<-pheno[,c(5,25)]
pheno<-cbind(pheno,do.call(rbind,stringr::str_split(pheno$tumour_stage,pattern = "T|N|M"))[,-1])


colnames(cl)
colnames(clinical)
colnames(pheno)

colnames(cl)<-c("ID","Gender","Status","Disease_status","Age")
colnames(clinical)<-c("ID","OS","OS.time")
colnames(pheno)<-c("ID","TNM","T","N","M")

cl<-cl[match(pheno$ID,cl$ID),]

cl$Age
cl$Age[nrow(cl)]
cl<-cl[-nrow(cl),] #-NA
pheno<-pheno[match(cl$ID,pheno$ID),]
clinical<-clinical[match(cl$ID,clinical$ID),]
cl<-cbind(cl,pheno[,-(1:2)],clinical[,-1])

exp<-read.csv("ICGC_lasso_cox.csv",header = T,row.names = 1,check.names = F)
cl<-cl[cl$ID%in%rownames(exp),]
sum(duplicated(cl$ID))
exp<-exp[rownames(exp)%in%cl$ID,]
cl<-cl[match(rownames(exp),cl$ID),]

rt<-cbind(exp[,c(1,2,(ncol(exp)-1),ncol(exp))],cl[,-c(1,3,4,9,10)])
rt<-rt[,-4]
colnames(rt)[3:ncol(rt)]<-c("Risk score","Gender","Age","T_stage","N_stage","M_stage")

table(rt$T_stage)
rt<-rt[!rt$T_stage=="X",]
rt$T_stage<-as.numeric(rt$T_stage)
table(rt$N_stage)
rt<-rt[!rt$N_stage=="X",]
rt$N_stage<-substr(rt$N_stage,1,1)
rt$N_stage<-as.numeric(rt$N_stage)

table(rt$M_stage)
rt<-rt[,-ncol(rt)]

table(rt$Gender)
rt$Gender<-ifelse(rt$Gender=="female",0,1)
rt$Gender<-as.numeric(rt$Gender)
rt$Age<-as.numeric(rt$Age)
uniTab=data.frame()
library(survival)
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(OStime, OSstatus) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="icgc.uniCox.txt",sep="\t",row.names=F,quote=F)


multiCox=coxph(Surv(OStime, OSstatus) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
rownames(multiTab)[1]<-"Risk score"
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="icgc.multiCox.txt",sep="\t",row.names=F,quote=F)

bioForest(coxFile="icgc.uniCox.txt",forestFile="icgc.uniForest.pdf",forestCol="#0D456B")
bioForest(coxFile="icgc.multiCox.txt",forestFile="icgc.multiForest.pdf",forestCol="#CC2216")
save(rt,file = "ICGC-riskscore-clincial.Rdata")

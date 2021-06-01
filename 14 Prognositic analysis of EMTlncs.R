# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
# TCGA
load("TCGA-riskscore-clincial.Rdata")
exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1,check.names = F)

exp<-exp[rownames(exp)%in%rownames(rt),]
rownames(exp)==rownames(rt)

rt<-cbind(exp,rt[,4:8])
uniTab=data.frame()
library(survival)
for(i in colnames(rt[,3:13])){
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
write.table(uniTab,file="tcga_model_lnc_uniCox.txt",sep="\t",row.names=F,quote=F)
bioForest=function(coxFile=null,forestFile=null){
  
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
  boxcolor = ifelse(as.numeric(hr) > 1, "#CC2216", "#0D456B")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
bioForest(coxFile="tcga_model_lnc_uniCox.txt",forestFile="tcga_model_lnc_uniForest.pdf")

# pheatmap
library(pheatmap)
rt<-rt[order(rt$riskScore,decreasing = F),]
annotation_col<-rt[,15:20]

annotation_col$Age<-ifelse(annotation_col$Age<60,"<60",">=60")
annotation_col$risk<-as.factor(annotation_col$risk)
annotation_col$Gender<-as.factor(annotation_col$Gender)
annotation_col$Age<-as.factor(annotation_col$Age)
annotation_col$Grade<-as.factor(annotation_col$Grade)
annotation_col$T_stage<-as.factor(annotation_col$T_stage)
annotation_col$N_stage<-as.factor(annotation_col$N_stage)

pheatmap(t(rt[,3:13]),cluster_rows = F,cluster_cols = F,show_colnames =  F,annotation_col=annotation_col)

# EMTlnc Prognostic analysis
library(survival) 
library(survminer) 
i<-colnames(rt)[12]
for (i in colnames(rt)[3:13]) {
  tiff(paste0(i,"_OS.tiff"),width = 480,height = 480)
  sur<-rt[,c("OStime","OSstatus",i)]
  colnames(sur)[3]<-"gene"
  sur$group<-ifelse(sur$gene>median(sur$gene),"High expression","Low expression")
  fit <- survfit(Surv(OStime, OSstatus)~group,data = sur)
  ggsurvplot(fit,sur,risk.table=TRUE,
           conf.int=TRUE,
           palette = c("#FF1E24","#0099FF"),
           pval=TRUE,
           pval.method=TRUE)
  dev.off()
}

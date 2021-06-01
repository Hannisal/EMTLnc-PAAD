# time:2021/4/7
rm(list = ls())
options(stringsAsFactors = F)
library("glmnet")
library("survival")
library("forestplot")
# TCGA
load("TCGA_exp_lnc.Rdata")
load("TCGA_clinical_Match.Rdata")

EMTlnc<-read.csv("TCGA_unicox_venn.csv",header =T)
exp_ATGlnc<-exp_lnc[EMTlnc$Genesymbol,]

OS<-clinical[,c(1,11,14)]
colnames(OS)<-c("ID","OSstatus","OStime")
OS$OSstatus<-ifelse(OS$OSstatus=="Dead",1,0)

OS<-OS[match(colnames(exp_ATGlnc),OS$ID),]
exp_ATGlnc<-t(exp_ATGlnc)
lasso_cox<-cbind(OS,exp_ATGlnc)

# lasso
x=as.matrix(lasso_cox[,c(4:ncol(lasso_cox))])
y=data.matrix(Surv(lasso_cox$OStime,lasso_cox$OSstatus))
set.seed(111)
fit=glmnet(x, y, family = "cox")

cvfit=cv.glmnet(x, y, family="cox",nfolds = 10)

plot(fit,xvar="lambda",label=T)
plot(cvfit)




coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)


trainFinalGeneExp=lasso_cox[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("ID","OStime","OSstatus",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(lasso_cox[,outCol],riskScore=as.vector(trainScore),risk)
write.csv(outTab,file = "TCGA_lasso_cox.csv",row.names = F)

# ICGC
load("ICGC_exp_lnc.Rdata")
load("clinical.Rdata")
rt<-t(exp_lnc)
rt=rt[,lassoGene]

OS<-clinical
colnames(OS)<-c("ID","OSstatus","OStime")
OS<-OS[!OS$OStime<30,]
rownames(rt)<-pheno$icgc_donor_id
rt<-rt[OS$ID,]
rt<-cbind(OS,rt)

testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("ID","OStime","OSstatus",lassoGene)
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.csv(outTab,file = "ICGC_lasso_cox.csv",row.names = F)
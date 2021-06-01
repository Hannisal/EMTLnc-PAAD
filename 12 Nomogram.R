# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
load("TCGA-riskscore-clincial.Rdata")

library(rms)
colnames(rt)[3]<-"Riskscore"
dd<-datadist(rt)
options(datadist="dd")
options(na.action="na.delete")
summary(rt$OStime)
coxpbc<-cph(formula = Surv(OStime,OSstatus) ~  Gender+Age+Grade+T_stage+N_stage+Riskscore ,data=rt,x=T,y=T,surv = T,na.action=na.delete)
print(coxpbc)
surv<-Survival(coxpbc) 
surv1<-function(x) surv(365,x)
surv2<-function(x) surv(730,x)
surv3<-function(x) surv(1095,x)

x<-nomogram(coxpbc,fun = list(surv1,surv2,surv3),lp=T,
            funlabel = c('1-year OS','2-year OS','3-year OS'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("TCGA-nomogram_classical.pdf",width = 16,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()
# ICGC
rm(list = ls())
load("ICGC-riskscore-clincial.Rdata")

library(rms)
colnames(rt)[3]<-"Riskscore"
dd<-datadist(rt)
options(datadist="dd")
options(na.action="na.delete")
summary(rt$OStime)
coxpbc<-cph(formula = Surv(OStime,OSstatus) ~  Gender+Age+T_stage+N_stage+Riskscore ,data=rt,x=T,y=T,surv = T,na.action=na.delete)
print(coxpbc)
surv<-Survival(coxpbc) 
surv1<-function(x) surv(365,x)
surv2<-function(x) surv(730,x)
surv3<-function(x) surv(1095,x)

x<-nomogram(coxpbc,fun = list(surv1,surv2,surv3),lp=T,
            funlabel = c('1-year OS','2-year OS','3-year OS'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("ICGC-nomogram_classical.pdf",width = 16,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

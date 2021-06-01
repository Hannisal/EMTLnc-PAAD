# time:2021/4/21
rm(list = ls())
options(stringsAsFactors = F)
# TCGA
load("TCGA-riskscore-clincial.Rdata")
library(survival)
library(rms)
colnames(rt)[3]<-"Riskscore"
for (i in 1:3) {
  pdf(paste0("TCGA Nomogram of ",i," Year.pdf"),height = 10,width = 10)
  res.cox1<-cph(Surv(OStime,OSstatus)~Riskscore,rt,surv = T,x=T,y=T,time.inc = i*365)
  cal<-calibrate(res.cox1,cmethod = "KM",method = "boot",u=i*365,m=46,B=138)
  plot(cal,lwd=2,lty=1,
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
       xlab = "Nomogram-Predicted Probability of OS",
       ylab = "Actual OS(proportion)",
       col=c(rgb(192,98,83,maxColorValue = 255)))
  dev.off()
}
# IGCG
load("ICGC-riskscore-clincial.Rdata")
library(survival)
library(rms)
colnames(rt)[3]<-"Riskscore"
for (i in 1:3) {
  pdf(paste0("ICGC Nomogram of ",i," Year.pdf"),height = 10,width = 10)
  res.cox1<-cph(Surv(OStime,OSstatus)~Riskscore,rt,surv = T,x=T,y=T,time.inc = i*365)
  cal<-calibrate(res.cox1,cmethod = "KM",method = "boot",u=i*365,m=28,B=89)
  plot(cal,lwd=2,lty=1,
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
       xlab = "Nomogram-Predicted Probability of OS",
       ylab = "Actual OS(proportion)",
       col=c(rgb(192,98,83,maxColorValue = 255)))
  dev.off()
}

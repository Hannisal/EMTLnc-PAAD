# time:2021/5/1
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
rt$Grade<-ifelse(rt$Grade=="G1","G1+G2",ifelse(rt$Grade=="G2","G1+G2","G3+G4"))
# G1+G2 G3+G4 
# 97    45 
table(rt$T_stage)
rt<-rt[!rt$T_stage=="X",]
rt$T_stage<-paste0("T",rt$T_stage)
rt$T_stage<-ifelse(rt$T_stage=="T1","T1+T2",ifelse(rt$T_stage=="T1","T1+T2","T3+T4"))
# T1+T2 T3+T4 
# 4   137 
table(rt$N_stage)
rt<-rt[!rt$N_stage=="X",]
rt$N_stage<-substr(rt$N_stage,1,1)
rt$N_stage<-paste0("N",rt$N_stage)
# N0  N1 
# 38 100 
table(rt$M_stage)
rt<-rt[,-ncol(rt)]

table(rt$Gender)
# FEMALE   MALE 
# 60     78 
rt$Age<-as.numeric(rt$Age)
rt$Age<-ifelse(rt$Age<60,"<60",">=60")
table(rt$Age)
# <60 >=60 
# 43   95
table(rt$Grade)
table(rt$T_stage)
table(rt$N_stage)
table(rt$M_stage)
table(rt$Gender)
table(rt$Age)
# G1+G2 G3+G4 
# 93    45 
# T1+T2 T3+T4 
# 4   134 
# N0  N1 
# 38 100
# FEMALE   MALE 
# 60     78
# <60 >=60 
# 43   95
colnames(rt)[3]<-"RiskScore"
# boxplot
library(ggpubr)
my_comparisons=list(group=levels(factor(rt$Gender)))
boxplot=ggboxplot(rt[,3:4], x="Gender", y="RiskScore", color="Gender",
                  xlab="Gender",
                  ylab="RiskScore",
                  legend.title="Gender",
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)
pdf(file=paste0("Gender",".","RiskScore",".pdf"),width=3.5,height=5)
print(boxplot)
dev.off()

# for
for (i in colnames(rt)[4:ncol(rt)]) {
  dat<-rt[,c("RiskScore",i)]
  my_comparisons=list(group=levels(factor(dat[,2])))
  boxplot=ggboxplot(dat, x=i, y="RiskScore", color=i,
                    xlab=i,
                    ylab="RiskScore",
                    legend.title=i,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  pdf(file=paste0(i,".","RiskScore",".pdf"),width=3,height=5)
  print(boxplot)
  dev.off()
}

# Stratification-survival analysis
library(survival) 
library(survminer) 
table(rt$Gender)
table(rt$Age)
table(rt$Grade)
table(rt$T_stage)
table(rt$N_stage)
# number of patients is >30,T1+T2 is excluded.
for (i in colnames(rt)[4:ncol(rt)]) {
  sur<-rt[,c("OStime","OSstatus","RiskScore",i)]
  for (j in unique(sur[,4])) {
    data<-sur[sur[,4]==j,]
    data$group<-ifelse(data$RiskScore>median(data$RiskScore),"High risk","Low risk")
    fit <- survfit(Surv(OStime, OSstatus)~group,data = data)
    surplot<-ggsurvplot(fit,data,risk.table=TRUE,
               conf.int=TRUE,
               palette = c("#FF1E24","#0099FF"),
               pval=TRUE,
               pval.method=TRUE)
    tiff(paste0(i,j,"_OS.tiff"),width = 480,height = 480)
    print(surplot)
    dev.off()
  }
}

# Age
i<-colnames(rt)[5]
sur<-rt[,c("OStime","OSstatus","RiskScore",i)]

j<-unique(sur[,4])[2]# change to 2 for <60
data<-sur[sur[,4]==j,]
data$group<-ifelse(data$RiskScore>median(data$RiskScore),"High risk","Low risk")
fit <- survfit(Surv(OStime, OSstatus)~group,data = data)
surplot<-ggsurvplot(fit,data,risk.table=TRUE,
                    conf.int=TRUE,
                    palette = c("#FF1E24","#0099FF"),
                    pval=TRUE,
                    pval.method=TRUE)
tiff(paste0(i,"lower60","_OS.tiff"),width = 480,height = 480)
print(surplot)
dev.off()

# ICGC
rm(list = ls())
load("ICGC-riskscore-clincial.Rdata")
exp<-read.csv("ICGC_lasso_cox.csv",header = T,row.names = 1,check.names = F)

exp<-exp[rownames(exp)%in%rownames(rt),]
rownames(exp)==rownames(rt)

rt<-cbind(exp,rt[,4:7])
rt$N_stage
rt<-rt[-62,] # remove na in T stage and N stage
table(rt$T_stage)
rt<-rt[!rt$T_stage=="X",]
rt$T_stage<-paste0("T",rt$T_stage)
rt$T_stage<-ifelse(rt$T_stage=="T1","T1+T2",ifelse(rt$T_stage=="T1","T1+T2","T3+T4"))
# T1+T2 T3+T4 
# 3    85
table(rt$N_stage)
rt<-rt[!rt$N_stage=="X",]
rt$N_stage<-substr(rt$N_stage,1,1)
rt$N_stage<-paste0("N",rt$N_stage)

table(rt$T_stage)
table(rt$N_stage)
table(rt$Gender)
table(annotation_col$Age)
# T1+T2 T3+T4 
# 3    85
# N0 N1 
# 30 58

# 0  1 
# 41 47 
# pheatmap
library(pheatmap)
rt<-rt[order(rt$riskScore,decreasing = F),]
annotation_col<-rt[,15:19]

annotation_col$Age<-ifelse(annotation_col$Age<60,"<60",">=60")
annotation_col$risk<-as.factor(annotation_col$risk)
annotation_col$Gender<-as.factor(annotation_col$Gender)
annotation_col$Age<-as.factor(annotation_col$Age)
annotation_col$T_stage<-as.factor(annotation_col$T_stage)
annotation_col$N_stage<-as.factor(annotation_col$N_stage)

pheatmap(t(rt[,3:13]),cluster_rows = F,cluster_cols = F,show_colnames =  F,annotation_col=annotation_col)

rt<-rt[,c(1:2,14:19)]
rt<-rt[,-4]
colnames(rt)[3]<-"RiskScore"
rt$Age<-ifelse(rt$Age<60,"<60",">=60")
# boxplot
library(ggpubr)
my_comparisons=list(group=levels(factor(rt$Gender)))
boxplot=ggboxplot(rt[,3:4], x="Gender", y="RiskScore", color="Gender",
                  xlab="Gender",
                  ylab="RiskScore",
                  legend.title="Gender",
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)
pdf(file=paste0("Gender",".","ICGC-RiskScore",".pdf"),width=3.5,height=5)
print(boxplot)
dev.off()

# for
for (i in colnames(rt)[4:ncol(rt)]) {
  dat<-rt[,c("RiskScore",i)]
  my_comparisons=list(group=levels(factor(dat[,2])))
  boxplot=ggboxplot(dat, x=i, y="RiskScore", color=i,
                    xlab=i,
                    ylab="RiskScore",
                    legend.title=i,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  pdf(file=paste0(i,".","ICGC-RiskScore",".pdf"),width=3,height=5)
  print(boxplot)
  dev.off()
}

# Stratification-survival analysis
library(survival) 
library(survminer) 
table(rt$Gender)
table(rt$Age)
table(rt$T_stage)
table(rt$N_stage)
# number of patients is >20,T1+T2 is excluded.
for (i in colnames(rt)[4:ncol(rt)]) {
  sur<-rt[,c("OStime","OSstatus","RiskScore",i)]
  for (j in unique(sur[,4])) {
    data<-sur[sur[,4]==j,]
    data$group<-ifelse(data$RiskScore>median(data$RiskScore),"High risk","Low risk")
    fit <- survfit(Surv(OStime, OSstatus)~group,data = data)
    surplot<-ggsurvplot(fit,data,risk.table=TRUE,
                        conf.int=TRUE,
                        palette = c("#FF1E24","#0099FF"),
                        pval=TRUE,
                        pval.method=TRUE)
    tiff(paste0(i,j,"ICGC_OS.tiff"),width = 480,height = 480)
    print(surplot)
    dev.off()
  }
}

# Age
i<-colnames(rt)[5]
sur<-rt[,c("OStime","OSstatus","RiskScore",i)]

j<-unique(sur[,4])[2]# change to 2 for <60
data<-sur[sur[,4]==j,]
data$group<-ifelse(data$RiskScore>median(data$RiskScore),"High risk","Low risk")
fit <- survfit(Surv(OStime, OSstatus)~group,data = data)
surplot<-ggsurvplot(fit,data,risk.table=TRUE,
                    conf.int=TRUE,
                    palette = c("#FF1E24","#0099FF"),
                    pval=TRUE,
                    pval.method=TRUE)
tiff(paste0(i,"lower60","ICGC_OS.tiff"),width = 480,height = 480)
print(surplot)
dev.off()

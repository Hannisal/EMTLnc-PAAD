# time:2021/3/26
rm(list = ls())
options(stringsAsFactors = F)
# Expression profiles
exp<-read.table("nature16965-s2_normalised_exp.txt",header = T)
# Clinical profiles
clinical<-read.table("donor.all_projects.overallSurvival_transfer_specimen",header = T,sep = "\t")
pheno<-read.table("specimen.all_projects.tsv",header = T,sep = "\t")
pheno<-pheno[pheno$submitted_donor_id%in%colnames(exp),]
sum(duplicated(pheno$submitted_donor_id))
pheno<-pheno[!duplicated(pheno$submitted_donor_id),]
exp<-exp[,pheno$submitted_donor_id]

clinical<-clinical[,-1];clinical<-clinical[!duplicated(clinical$X_PATIENT),]
clinical<-clinical[clinical$X_PATIENT%in%pheno$icgc_donor_id,]
# Annotation
load("anno.Rdata")
lnc_anno$gene_id<-substr(lnc_anno$gene_id,1,15)
mRNA_anno$gene_id<-substr(mRNA_anno$gene_id,1,15)
exp_lnc<-exp[rownames(exp)%in%lnc_anno$gene_id,]
exp_mRNA<-exp[rownames(exp)%in%mRNA_anno$gene_id,]

library(dplyr)
library(tibble)
id<-mRNA_anno[,1:2]
colnames(id)<-c("gene_name","gene_id")
exp_mRNA<-exp_mRNA %>%
  rownames_to_column(var = "gene_id")%>%
  inner_join(id,by = "gene_id")%>%
  select(-gene_id)%>%
  select(gene_name,everything())%>%
  mutate(rowMean = rowMeans(.[grep("ENSG",names(.))]))%>%
  arrange(desc(rowMean))%>%
  distinct(gene_name,.keep_all = T)%>%
  select(-rowMean)%>%
  column_to_rownames(var = "gene_name")

library(dplyr)
library(tibble)
id<-lnc_anno[,1:2]
colnames(id)<-c("gene_name","gene_id")
exp_lnc<-exp_lnc %>%
  rownames_to_column(var = "gene_id")%>%
  inner_join(id,by = "gene_id")%>%
  select(-gene_id)%>%
  select(gene_name,everything())%>%
  mutate(rowMean = rowMeans(.[grep("ENSG",names(.))]))%>%
  arrange(desc(rowMean))%>%
  distinct(gene_name,.keep_all = T)%>%
  select(-rowMean)%>%
  column_to_rownames(var = "gene_name")
nrow(exp_lnc)
# 1440
ncol(exp_lnc)
# 95
nrow(exp_mRNA)
# 15146
save(exp_mRNA,file = "ICGC_exp_mRNA.Rdata")
save(exp_lnc,file = "ICGC_exp_lnc.Rdata")
save(pheno,clinical,file = "clinical.Rdata")

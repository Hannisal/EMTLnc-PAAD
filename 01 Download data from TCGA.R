# time:2021/2/20
rm(list = ls())
options(stringsAsFactors = F)
cancer_typer="TCGA-PAAD"
# Expression profiles
x = read.table("expdata/0ae4ff1f-e2d3-46e0-95a2-0ea80a4ebb63/574df2fc-a608-49c5-8e83-f26d03ef8bb3.FPKM.txt.gz")
x2 = read.table("expdata/0c2840a2-3a49-4f22-ae21-1cfbb0034212/fef65b57-c58d-4050-8de4-f09f5cd616ce.FPKM.txt.gz")
identical(x$V1,x2$V1)

count_files = dir("expdata/",pattern = "*.FPKM.txt.gz$",recursive = T)

ex = function(x){
  result <- read.table(file.path("expdata/",x),row.names = 1,sep = "\t") 
  return(result)
}
exp = lapply(count_files,ex)
exp <- do.call(cbind,exp)
dim(exp)
exp[1:4,1:4]

meta <- jsonlite::fromJSON("metadata.cart.exp.json")
colnames(meta)
meta$associated_entities[1]

ids <- meta$associated_entities;class(ids)
ids[[1]]

class(ids[[1]][,1])


colnames(ids)
ID = sapply(ids,function(x){x[,1]})
file2id = data.frame(file_name = meta$file_name, ID = ID)


count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
count_files2[1] %in% file2id$file_name
file2id = file2id[match(count_files2,file2id$file_name),]
colnames(exp) = file2id$ID
exp[1:4,1:4]
exp<-log2(exp+1)
save(exp,file = "exp_TCGA_FPKM.Rdata")

# Clinical profiles
library(XML)
result <- xmlParse("./clinical/0cffa300-8b92-4424-97a4-32ec9584e95f/nationwidechildrens.org_clinical.TCGA-2J-AABH.xml")
rootnode <- xmlRoot(result)
rootsize <- xmlSize(rootnode)
print(rootnode[1])

xmldataframe <- xmlToDataFrame(rootnode[2])
head(t(xmlToDataFrame(rootnode[2])))

xmls = dir("clinical/",pattern = "*.xml$",recursive = T)

td = function(x){
  result <- xmlParse(file.path("clinical/",x))
  rootnode <- xmlRoot(result)
  xmldataframe <- xmlToDataFrame(rootnode[2])
  return(t(xmldataframe))
}

cl = lapply(xmls,td)
cl_df <- t(do.call(cbind,cl))
cl_df[1:3,1:3]

clinical = data.frame(cl_df)
clinical[1:4,1:4]

save(clinical,file = "clinical_TCGA.Rdata")


load("clinical_TCGA.Rdata")
clinical<-clinical[,c(4,12:15,19,21,26,35,39:42)]
head(clinical)

clinical$OStime<-ifelse(clinical$days_to_last_followup=="",clinical$days_to_death,clinical$days_to_last_followup)
## OS missing and OS<30
clinical$OStime<-as.numeric(clinical$OStime)
clinical<-clinical[!clinical$OStime=="",]
clinical<-clinical[!clinical$OStime<30,]

load("exp_TCGA_FPKM.Rdata")


library(stringr)
table(str_sub(colnames(exp),14,15))

group_list = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)
## tumor
exp<-exp[,group_list=="tumor"]
sum(duplicated(substr(colnames(exp),1,12)))

a<-colnames(exp)[duplicated(substr(colnames(exp),1,12))]
colnames(exp)[substr(colnames(exp),1,12)=="TCGA-HZ-A9TJ"]
## Choose:"TCGA-HZ-A9TJ-01A-11R-A41I-07"

exp<-exp[,!colnames(exp)=="TCGA-HZ-A9TJ-06A-11R-A41B-07"]
colnames(exp)<-substr(colnames(exp),1,12)
exp<-exp[,colnames(exp)%in%clinical$bcr_patient_barcode]
clinical<-clinical[clinical$bcr_patient_barcode%in%colnames(exp),]

clinical<-clinical[match(colnames(exp),clinical$bcr_patient_barcode),]

save(exp,file = "TCGA_FPKM_Match.Rdata")
save(clinical,file = "TCGA_clinical_Match.Rdata")

# Annotation
options(stringsAsFactors = F)
if(!file.exists("anno.Rdata")){
  #BiocManager::install("rtracklayer")
  library(rtracklayer)
  x = rtracklayer::import("gencode.v22.annotation.gtf")
  class(x)
  x2 = as.data.frame(x);dim(x2)
  colnames(x2)
  tj = as.data.frame(table(x2$type));tj
  
  tj2 = as.data.frame(table(x2$gene_type));tj2
  
  nrow(x2);x2 = x2[x2$type=="gene",];nrow(x2)
  tj3 = as.data.frame(table(x2$gene_type));tj3
  
  lnc_bype = c("3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", "non_coding", "processed_transcript", "sense_intronic" , "sense_overlapping")
  table(x2$gene_type %in% lnc_bype)
  table(x2$gene_type == "protein_coding")
  lnc_anno = x2[x2$gene_type %in% lnc_bype,c("gene_name","gene_id","gene_type")]
  mRNA_anno = x2[x2$gene_type == "protein_coding",c("gene_name","gene_id","gene_type")]
  save(lnc_anno,mRNA_anno,file = "anno.Rdata")
}
load("anno.Rdata")


load("TCGA_FPKM_Match.Rdata")
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
# 14805
ncol(exp_lnc)
# 146
nrow(exp_mRNA)
# 19712
save(exp_mRNA,file = "TCGA_exp_mRNA.Rdata")
save(exp_lnc,file = "TCGA_exp_lnc.Rdata")

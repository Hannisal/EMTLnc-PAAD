# time:2021/5/1
rm(list = ls())
options(stringsAsFactors = F)
# TCGA
load("TCGA_exp_mRNA.Rdata")
exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1)
group_list<-as.factor(exp$risk)
group_list<-factor(group_list,levels = c("low","high") )


table(group_list)
fdrFilter=0.05
logFCfilter=1


outTab=data.frame()
dimnames=list(rownames(exp_mRNA),colnames(exp_mRNA))
data=matrix(as.numeric(as.matrix(exp_mRNA)),nrow=nrow(exp_mRNA),dimnames=dimnames)
data=2^data-1

for(i in row.names(data)){
  geneName=i
  rt=rbind(expression=data[i,],group_list=group_list)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ group_list, data=rt)
  conGeneMeans=mean(data[i,group_list=="low"])
  treatGeneMeans=mean(data[i,group_list=="high"])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,group_list=="low"])
  treatMed=median(data[i,group_list=="high"])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
outDiff_gene<-outDiff$gene
write.table(outDiff_gene,file = "input_gene.txT",quote = F,row.names = F)

rm(list = ls())
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
# prepare for GSEA
load("TCGA_exp_mRNA.Rdata")
exp<-read.csv("TCGA_lasso_cox.csv",header = T,row.names = 1)
group_list<-as.factor(exp$risk)
levels(group_list)<-c("low","high")

gmtFile="h.all.v7.4.symbols.gmt"
gmt=read.gmt(gmtFile)

dimnames=list(rownames(exp_mRNA),colnames(exp_mRNA))
data=matrix(as.numeric(as.matrix(exp_mRNA)),nrow=nrow(exp_mRNA),dimnames=dimnames)
data=2^data-1

col1<-data.frame(NAME=rownames(data),Description="mRNA")

input<-cbind(col1,data)

write.table(input,file = "GSEA_input.txt",sep = "\t",row.names = F,quote = F)

paste0(exp$risk,collapse = " ")
write.csv(paste0(exp$risk,collapse = " "),file = "group_GSEA.csv")

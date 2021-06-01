# time:2021/4/7
rm(list = ls())
options(stringsAsFactors = F)
load("TCGA_exp_lnc.Rdata")
gsva<-read.csv("gsva_output.csv",header = T,row.names = 1,check.names = F)
rownames(gsva)
HallMarks<-c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gsva<-gsva[HallMarks,]

mydata<-exp_lnc
datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0

library(WGCNA)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
## cutoff:70
clust = cutreeStatic(sampleTree, cutHeight = 750, minSize = 10)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
# gene module
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("1Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
# gene module
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

samples=t(gsva)
samples<-samples[rownames(samples)%in%rownames(datExpr),]

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)

modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
rownames(textMatrix)<-colnames(MEsWW)
colnames(textMatrix)<-"Cor_pvalue"
write.csv(textMatrix,file = "WGCNA model cor and p.csv")
dim(textMatrix) = dim(modTraitCor)
# abs(cor)>0.5&p.value<0.05
# MElightcyan
# MEturquoise 
# names (colors) of the modules
modNames = substring(names(MEsWW), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEsWW, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


geneTraitSignificance = as.data.frame(cor(datExpr, samples, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", "EMT", sep="")
names(GSPvalue) = paste("p.GS.", "EMT", sep="")

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EMT",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "#0099FF",pch = 16)

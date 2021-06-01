# time:2021/4/7
rm(list = ls())
options(stringsAsFactors = F)
library(msigdbr)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)
load("TCGA_exp_mRNA.Rdata")

msigdbr_show_species()
h <- msigdbr(species = "Homo sapiens",
             category = "H")


h <- select(h, gs_name, gene_symbol) %>%
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol))


save(h, file = "hallmark.gs.RData")

gsva_es <- gsva(as.matrix(exp_mRNA), h)

write.csv(gsva_es, "gsva_output.csv", quote = F)

rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

if(T){
  library(readr)
  library(readxl)
  library(tidyverse)
  library(GEOquery)
  library(plyr)
  library(circlize)
  library(ComplexHeatmap)
  library(devtools)
  library(SummarizedExperiment)
  library(parallel)
  library(KEGG.db)
  library(data.table)
  library(edgeR)
  library(limma)
  library(TCGAbiolinks)
  library(GEOquery)
  library(readxl)
  library(clusterProfiler)
  library(openxlsx)
  library(org.Hs.eg.db)
}

if(!dir.exists('02_PathwayEnrichment')){
  dir.create('02_PathwayEnrichment')
}

allDiff<-read.csv(file = '01_DEGs/SAE vs Ctrl_DEGs.csv',check.names = F)
colnames(allDiff)[1]<-"SYMBOL"
ID <- bitr(allDiff$SYMBOL, fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db) 
ID<-as.data.frame(ID)
id<-intersect(ID$SYMBOL,allDiff$SYMBOL)
ID<-ID[ID$SYMBOL%in%id,]
ID<-ID[order(ID$SYMBOL,ID$ENTREZID,decreasing = T),]
ID<-ID[!duplicated(ID$SYMBOL),]
allDiff<-allDiff[allDiff$SYMBOL%in%id,]
df<-left_join(ID,allDiff, by = c("SYMBOL"))

sig<-df[df$logFC>0&df$P.Value<0.05,]
gene<-sig$ENTREZID

enrich <- enrichKEGG(gene=gene, organism='hsa', keyType="kegg", use_internal_data=T)
head(enrich)
kegg<-as.data.frame(enrich@result)
write.xlsx(kegg,file = '02_PathwayEnrichment/enrichKEGG_SAE vs Ctrl.xlsx')


##############
allDiff<-read.csv(file = '01_DEGs/SAE vs SP_DEGs.csv',check.names = F)
colnames(allDiff)[1]<-"SYMBOL"
ID <- bitr(allDiff$SYMBOL, fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db) 
ID<-as.data.frame(ID)
id<-intersect(ID$SYMBOL,allDiff$SYMBOL)
ID<-ID[ID$SYMBOL%in%id,]
ID<-ID[order(ID$SYMBOL,ID$ENTREZID,decreasing = T),]
ID<-ID[!duplicated(ID$SYMBOL),]
allDiff<-allDiff[allDiff$SYMBOL%in%id,]
df<-left_join(ID,allDiff, by = c("SYMBOL"))

sig<-df[df$logFC>0&df$P.Value<0.05,]

gene<-sig$ENTREZID
enrich <- enrichKEGG(gene=gene, organism='hsa', keyType="kegg", use_internal_data=T)
head(enrich)
kegg<-as.data.frame(enrich@result)
write.xlsx(kegg,file = '02_PathwayEnrichment/enrichKEGG_SAE vs SP.xlsx')



##############
allDiff<-read.csv(file = '01_DEGs/SP vs Ctrl_DEGs.csv',check.names = F)
colnames(allDiff)[1]<-"SYMBOL"
ID <- bitr(allDiff$SYMBOL, fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db) 
ID<-as.data.frame(ID)
id<-intersect(ID$SYMBOL,allDiff$SYMBOL)
ID<-ID[ID$SYMBOL%in%id,]
ID<-ID[order(ID$SYMBOL,ID$ENTREZID,decreasing = T),]
ID<-ID[!duplicated(ID$SYMBOL),]
allDiff<-allDiff[allDiff$SYMBOL%in%id,]
df<-left_join(ID,allDiff, by = c("SYMBOL"))

sig<-df[df$logFC>0&df$P.Value<0.05,]
gene<-sig$ENTREZID

enrich <- enrichKEGG(gene=gene, organism='hsa', keyType="kegg", use_internal_data=T)
head(enrich)
kegg<-as.data.frame(enrich@result)
write.xlsx(kegg,file = '02_PathwayEnrichment/enrichKEGG_SP vs Ctrl.xlsx')

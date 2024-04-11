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
  library(data.table)
  library(edgeR)
  library(limma)
  library(TCGAbiolinks)
  library(GEOquery)
  library(readxl)
  library(openxlsx)
}

if(!dir.exists('01_DEGs')){
  dir.create('01_DEGs')
}

df<-read_xlsx('YuanSF_RNASeq_CPOS-230919-SY-19861a_Gene_DiffExpressed.xlsx',sheet = 1)
df<-as.data.frame(df)
colnames(df)
exp<-df[,c('gene_name',colnames(df)[8:13])]
# rownames(exp)<-exp$gene_name
colnames(exp)
rownames(exp)<-paste0('a',1:nrow(exp))
dd<-data.frame(gene = exp$gene_name, index = paste0('a',1:nrow(exp)), num = c(0))
dd$num<-apply(exp[,-1],1,mean)
dd<-dd[order(dd$gene,dd$num,decreasing = T),]
dd1<-dd[!duplicated(dd$gene),]
exp<-exp[dd1$index,]
rownames(exp)<-exp$gene_name
exp<-exp[,-1]
colnames(exp)<-gsub(' TPM','',colnames(exp))
colnames(exp)
mut<-data.frame(row.names = colnames(exp),group= c(rep('other',3),rep('me',3)))
mut$group1<-mut$group
both<-intersect(colnames(exp),rownames(mut))
exp<-exp[,both]
cli<-mut[both,]

# exp<-log2(exp+1)
save(exp,cli,file = '01_DEGs/SP vs Ctrl.Rdata')

group_list= cli$group ##High vs Low

do_limma_array <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  fit <- lmFit(exprSet, design)
  group_list
  cont.matrix=makeContrasts(contrasts=c('me-other'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='me-other', n=Inf)
  DEG_limma = na.omit(tempOutput)
  head(DEG_limma)
  return(DEG_limma)
}


deg1=do_limma_array(exp,group_list)
# deg1=do_edgeR_counts(exp,group_list)
head(deg1)[,1:5]
deg1$GeneSymbol<-rownames(deg1)
deg1<-deg1[,c(7,1:6)]
colnames(deg1)
# colnames(deg1)[5:6]<-c("PValue","FDR")

deg1$group<-c('SP vs Ctrl')
write.csv(deg1,file = '01_DEGs/SP vs Ctrl_DEGs.csv',row.names = F,quote = F)


# rm(list = ls())
df<-read_xlsx('YuanSF_RNASeq_CPOS-230919-SY-19861a_Gene_DiffExpressed.xlsx',sheet = 2)
df<-as.data.frame(df)
colnames(df)
exp<-df[,c('gene_name',colnames(df)[8:13])]
# rownames(exp)<-exp$gene_name
rownames(exp)<-paste0('a',1:nrow(exp))
dd<-data.frame(gene = exp$gene_name, index = paste0('a',1:nrow(exp)), num = c(0))
dd$num<-apply(exp[,-1],1,mean)
dd<-dd[order(dd$gene,dd$num,decreasing = T),]
dd1<-dd[!duplicated(dd$gene),]
exp<-exp[dd1$index,]
rownames(exp)<-exp$gene_name
exp<-exp[,-1]
colnames(exp)<-gsub(' TPM','',colnames(exp))

mut<-data.frame(row.names = colnames(exp),group= c(rep('other',3),rep('me',3)))
mut$group1<-mut$group
both<-intersect(colnames(exp),rownames(mut))
exp<-exp[,both]
cli<-mut[both,]

# exp<-log2(exp+1)
save(exp,cli,file = '01_DEGs/SAE vs Ctrl.Rdata')

group_list= cli$group 

deg1=do_limma_array(exp,group_list)
# deg1=do_edgeR_counts(exp,group_list)
head(deg1)[,1:5]
deg1$GeneSymbol<-rownames(deg1)
deg1<-deg1[,c(7,1:6)]
colnames(deg1)
# colnames(deg1)[5:6]<-c("PValue","FDR")

deg1$group<-c('SAE vs Ctrl')
write.csv(deg1,file = '01_DEGs/SAE vs Ctrl_DEGs.csv',row.names = F,quote = F)


df<-read_xlsx('YuanSF_RNASeq_CPOS-230919-SY-19861a_Gene_DiffExpressed.xlsx',sheet = 3)
df<-as.data.frame(df)
colnames(df)
exp<-df[,c('gene_name',colnames(df)[8:13])]
# rownames(exp)<-exp$gene_name
rownames(exp)<-paste0('a',1:nrow(exp))
dd<-data.frame(gene = exp$gene_name, index = paste0('a',1:nrow(exp)), num = c(0))
dd$num<-apply(exp[,-1],1,mean)
dd<-dd[order(dd$gene,dd$num,decreasing = T),]
dd1<-dd[!duplicated(dd$gene),]
exp<-exp[dd1$index,]
rownames(exp)<-exp$gene_name
exp<-exp[,-1]
colnames(exp)<-gsub(' TPM','',colnames(exp))
colnames(exp)
mut<-data.frame(row.names = colnames(exp),group= c(rep('other',3),rep('me',3)))
mut$group1<-mut$group
both<-intersect(colnames(exp),rownames(mut))
exp<-exp[,both]
cli<-mut[both,]


save(exp,cli,file = '01_DEGs/SAE vs SP.Rdata')

group_list= cli$group 

deg1=do_limma_array(exp,group_list)
# deg1=do_edgeR_counts(exp,group_list)
head(deg1)[,1:5]
deg1$GeneSymbol<-rownames(deg1)
deg1<-deg1[,c(7,1:6)]
colnames(deg1)
# colnames(deg1)[5:6]<-c("PValue","FDR")

deg1$group<-c('SAE vs SP')
write.csv(deg1,file = '01_DEGs/SAE vs SP_DEGs.csv',row.names = F,quote = F)

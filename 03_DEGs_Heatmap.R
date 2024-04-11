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
  library(clusterProfiler)
  library(openxlsx)
  library(org.Hs.eg.db)
}

if(!dir.exists('03_DEGs_Heatmap')){
  dir.create('03_DEGs_Heatmap')
}

index<-read.xlsx('02_PathwayEnrichment/enrichKEGG_SAE vs Ctrl.xlsx')
index<-as.data.frame(index)

index$Description[grep('Apoptosis',index$Description)]
index$Description[grep('Cellular senescence',index$Description)]
index$Description[grep('Autophagy',index$Description)]
want<-c('Apoptosis','Apoptosis - multiple species"','Cellular senescence','Autophagy - animal','Autophagy - other')
ind<-index[index$Description%in%want,]
keggRes<-ind
keggRes$RepGenes <- NA_character_ 

pl<-list()
for (i in 1:nrow(keggRes)) {
  intGenes <- unlist(strsplit(keggRes$geneID[i], "/"))
  dd<-data.frame(path = keggRes$Description[i], geneid = intGenes)
  pl[[i]]<-dd
}

pl_df<-do.call(rbind,pl)
pl_df<-as.data.frame(pl_df)

ID <- bitr(pl_df$geneid, fromType = "ENTREZID",
           toType = c("SYMBOL"),
           OrgDb = org.Hs.eg.db) 
ID<-as.data.frame(ID)
class(pl_df$geneid);class(ID$ENTREZID)
both_id<-intersect(pl_df$geneid,ID$ENTREZID)
ID<-ID[ID$ENTREZID%in%both_id,]
pl_df<-pl_df[pl_df$geneid%in%both_id,]
colnames(pl_df)[2]<-c('ENTREZID')
pl_df<-left_join(pl_df,ID, by = c('ENTREZID'))

deg<-read.csv('01_DEGs/SAE vs Ctrl_DEGs.csv',check.names = F)
deg<-deg[deg$GeneSymbol%in%pl_df$SYMBOL,]
all(deg$P.Value<0.05)
rownames(deg)<-deg$GeneSymbol
deg1<-deg[ID$SYMBOL,]
colnames(deg1)
deg1<-deg1[,c('GeneSymbol','P.Value')]

deg<-read.csv('01_DEGs/SAE vs SP_DEGs.csv',check.names = F)
deg<-deg[deg$GeneSymbol%in%pl_df$SYMBOL,]
all(deg$P.Value<0.05)
rownames(deg)<-deg$GeneSymbol
deg2<-deg[ID$SYMBOL,]
colnames(deg2)
deg2<-deg2[,c('GeneSymbol','P.Value')]


deg<-read.csv('01_DEGs/SP vs Ctrl_DEGs.csv',check.names = F)
deg<-deg[deg$GeneSymbol%in%pl_df$SYMBOL,]
all(deg$P.Value<0.05)
rownames(deg)<-deg$GeneSymbol
deg3<-deg[ID$SYMBOL,]
colnames(deg3)
deg3<-deg3[,c('GeneSymbol','P.Value')]

colnames(deg1)[2]<-c('SAE_Ctrl')
colnames(deg2)[2]<-c('SAE_SP')
colnames(deg3)[2]<-c('SP_Ctrl')
deg<-left_join(deg1,deg2,by = c('GeneSymbol'))
deg<-left_join(deg,deg3,by = c('GeneSymbol'))
write.csv(deg,file = '03_DEGs_Heatmap/DEGs_Pvalue.csv',row.names = F,quote = F)

load('01_DEGs/SAE vs Ctrl.Rdata')
exp1<-exp[deg$GeneSymbol,]
load('01_DEGs/SAE vs SP.Rdata')
exp2<-exp[deg$GeneSymbol,]

colnames(exp1)
colnames(exp2)
exp<-cbind(exp1,exp2[,1:3])
colnames(exp)
exp<-exp[,c(1:3,7:9,4:6)]
write.csv(exp,file = '03_DEGs_Heatmap/exp.csv',row.names = T,quote = F)

path <- unique(pl_df$path)[1]
pl_df1<-pl_df[pl_df$path%in%path,]

exp<-log2(exp+1)

rownames(deg)<-deg$GeneSymbol
deg$SAE_Ctrl<-ifelse(deg$SAE_Ctrl<0.0001,'****',
                     ifelse(deg$SAE_Ctrl<0.001,'***',
                            ifelse(deg$SAE_Ctrl<0.01,'**',
                                   ifelse(deg$SAE_Ctrl<0.05,'*','ns'))))

deg$SAE_SP<-ifelse(deg$SAE_SP<0.0001,'****',
                   ifelse(deg$SAE_SP<0.001,'***',
                          ifelse(deg$SAE_SP<0.01,'**',
                                 ifelse(deg$SAE_SP<0.05,'*','ns'))))


deg$SP_Ctrl<-ifelse(deg$SP_Ctrl<0.0001,'****',
                    ifelse(deg$SP_Ctrl<0.001,'***',
                           ifelse(deg$SP_Ctrl<0.01,'**',
                                  ifelse(deg$SP_Ctrl<0.05,'*','ns'))))

ha<-deg[deg$GeneSymbol%in%pl_df1$SYMBOL,]
colnames(ha)
unique(ha$SAE_Ctrl)
c1<-c('#F2E5E3','#D2A49D','#B26357','#7E433A')
names(c1)<-c('*','**','***','****')
unique(ha$SAE_SP)

c2<-c('#F2E5E3','#D2A49D','#B26357','#7E433A','gray')
names(c2)<-c('*','**','***','****','ns')

unique(ha$SP_Ctrl)
c3<-c('#F2E5E3','gray')
names(c3)<-c('*','ns')
ha1 = rowAnnotation("SAE vs Ctrl"=ha$SAE_Ctrl,
                        "SAE vs SP"=ha$SAE_SP,
                        "SP vs Ctrl"=ha$SP_Ctrl,
                       col=list(
                         "SAE vs Ctrl"=c1,
                         "SAE vs SP"=c2,
                         "SP vs Ctrl"=c3
                       ),
                       show_annotation_name = TRUE,
                       # annotation_name_side="right",
                       annotation_name_gp = gpar(fontsize = 7),
                       show_legend=T,
                       simple_anno_size = unit(0.2, "cm"),
                       annotation_height=unit.c(unit(0.07, "cm")))

expp<-exp[pl_df1$SYMBOL,]
expp<-scale(as.matrix(t(expp)))
expp[expp>2]=2
expp[expp<c(-2)]<-c(-2)

expp<-as.matrix(t(expp))

library(circlize)
heat_colors2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = paste0('03_DEGs_Heatmap/',path,'.pdf'),width = 8,height = 12)
anti <- Heatmap(as.matrix(expp),
                name='Z-score',
                # top_annotation = ha,
                # cluster_rows = T,
                left_annotation = ha1,
                col=heat_colors2,
                color_space = "RGB",
                cluster_columns = FALSE,
                row_order=NULL,
                column_order= colnames(expp),
                show_column_names = T,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                gap = unit(1, "mm"),
                column_title_gp = gpar(fontsize = 7),
                width=unit(2.6, "cm"),
                height=unit(nrow(expp)*0.26, "cm"),
                show_heatmap_legend = T,
                heatmap_legend_param=list(labels_gp = gpar(fontsize = 7, fontface = "bold"),
                                          title_gp = gpar(fontsize = 7, fontface = "bold"))) 

anti
dev.off()

path <- unique(pl_df$path)[2]
pl_df1<-pl_df[pl_df$path%in%path,]

ha<-deg[deg$GeneSymbol%in%pl_df1$SYMBOL,]
colnames(ha)
unique(ha$SAE_Ctrl)
c1<-c('#F2E5E3','#D2A49D','#B26357','#7E433A')
names(c1)<-c('*','**','***','****')
unique(ha$SAE_SP)
c2<-c('#F2E5E3','#D2A49D','#B26357','#7E433A','gray')
names(c2)<-c('*','**','***','****','ns')

unique(ha$SP_Ctrl)
c3<-c('#F2E5E3','#D2A49D','#B26357','gray')
names(c3)<-c('*','**','***','ns')

ha1 = rowAnnotation("SAE vs Ctrl"=ha$SAE_Ctrl,
                    "SAE vs SP"=ha$SAE_SP,
                    "SP vs Ctrl"=ha$SP_Ctrl,
                    col=list(
                      "SAE vs Ctrl"=c1,
                      "SAE vs SP"=c2,
                      "SP vs Ctrl"=c3
                    ),
                    show_annotation_name = TRUE,
                    # annotation_name_side="right",
                    annotation_name_gp = gpar(fontsize = 7),
                    show_legend=T,
                    simple_anno_size = unit(0.2, "cm"),
                    annotation_height=unit.c(unit(0.07, "cm")))

expp<-exp[pl_df1$SYMBOL,]
expp<-scale(as.matrix(t(expp)))
expp[expp>2]=2
expp[expp<c(-2)]<-c(-2)

expp<-as.matrix(t(expp))

library(circlize)
heat_colors2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = paste0('03_DEGs_Heatmap/',path,'.pdf'),width = 8,height = 12)
anti <- Heatmap(as.matrix(expp),
                name='Z-score',
                # top_annotation = ha,
                # cluster_rows = T,
                left_annotation = ha1,
                col=heat_colors2,
                color_space = "RGB",
                cluster_columns = FALSE,
                row_order=NULL,
                column_order= colnames(expp),
                show_column_names = T,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                gap = unit(1, "mm"),
                column_title_gp = gpar(fontsize = 7),
                width=unit(2.6, "cm"),
                height=unit(nrow(expp)*0.26, "cm"),
                show_heatmap_legend = T,
                heatmap_legend_param=list(labels_gp = gpar(fontsize = 7, fontface = "bold"),
                                          title_gp = gpar(fontsize = 7, fontface = "bold"))) 

anti
dev.off()


path <- unique(pl_df$path)[3]
pl_df1<-pl_df[pl_df$path%in%path,]

ha<-deg[deg$GeneSymbol%in%pl_df1$SYMBOL,]
colnames(ha)
unique(ha$SAE_Ctrl)
c1<-c('#F2E5E3','#D2A49D','#B26357','#7E433A')
names(c1)<-c('*','**','***','****')
unique(ha$SAE_SP)
c2<-c('#F2E5E3','#D2A49D','#B26357','#7E433A','gray')
names(c2)<-c('*','**','***','****','ns')

unique(ha$SP_Ctrl)
c3<-c('#F2E5E3','#D2A49D','#B26357','gray')
names(c3)<-c('*','**','***','ns')

ha1 = rowAnnotation("SAE vs Ctrl"=ha$SAE_Ctrl,
                    "SAE vs SP"=ha$SAE_SP,
                    "SP vs Ctrl"=ha$SP_Ctrl,
                    col=list(
                      "SAE vs Ctrl"=c1,
                      "SAE vs SP"=c2,
                      "SP vs Ctrl"=c3
                    ),
                    show_annotation_name = TRUE,
                    # annotation_name_side="right",
                    annotation_name_gp = gpar(fontsize = 7),
                    show_legend=T,
                    simple_anno_size = unit(0.2, "cm"),
                    annotation_height=unit.c(unit(0.07, "cm")))

expp<-exp[pl_df1$SYMBOL,]
expp<-scale(as.matrix(t(expp)))
expp[expp>2]=2
expp[expp<c(-2)]<-c(-2)

expp<-as.matrix(t(expp))

library(circlize)
heat_colors2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = paste0('03_DEGs_Heatmap/',path,'.pdf'),width = 8,height = 12)
anti <- Heatmap(as.matrix(expp),
                name='Z-score',
                # top_annotation = ha,
                # cluster_rows = T,
                left_annotation = ha1,
                col=heat_colors2,
                color_space = "RGB",
                cluster_columns = FALSE,
                row_order=NULL,
                column_order= colnames(expp),
                show_column_names = T,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                gap = unit(1, "mm"),
                column_title_gp = gpar(fontsize = 7),
                width=unit(2.6, "cm"),
                height=unit(nrow(expp)*0.26, "cm"),
                show_heatmap_legend = T,
                heatmap_legend_param=list(labels_gp = gpar(fontsize = 7, fontface = "bold"),
                                          title_gp = gpar(fontsize = 7, fontface = "bold"))) 

anti
dev.off()


path <- unique(pl_df$path)[4]
pl_df1<-pl_df[pl_df$path%in%path,]

ha<-deg[deg$GeneSymbol%in%pl_df1$SYMBOL,]
colnames(ha)
unique(ha$SAE_Ctrl)
c1<-c('#F2E5E3','#D2A49D','#B26357')
names(c1)<-c('*','**','***')
unique(ha$SAE_SP)
c2<-c('#F2E5E3','#D2A49D','#B26357','gray')
names(c2)<-c('*','**','***','ns')

unique(ha$SP_Ctrl)
c3<-c('gray')
names(c3)<-c('ns')

ha1 = rowAnnotation("SAE vs Ctrl"=ha$SAE_Ctrl,
                    "SAE vs SP"=ha$SAE_SP,
                    "SP vs Ctrl"=ha$SP_Ctrl,
                    col=list(
                      "SAE vs Ctrl"=c1,
                      "SAE vs SP"=c2,
                      "SP vs Ctrl"=c3
                    ),
                    show_annotation_name = TRUE,
                    # annotation_name_side="right",
                    annotation_name_gp = gpar(fontsize = 7),
                    show_legend=T,
                    simple_anno_size = unit(0.2, "cm"),
                    annotation_height=unit.c(unit(0.07, "cm")))

expp<-exp[pl_df1$SYMBOL,]
expp<-scale(as.matrix(t(expp)))
expp[expp>2]=2
expp[expp<c(-2)]<-c(-2)

expp<-as.matrix(t(expp))

library(circlize)
heat_colors2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = paste0('03_DEGs_Heatmap/',path,'.pdf'),width = 8,height = 12)
anti <- Heatmap(as.matrix(expp),
                name='Z-score',
                # top_annotation = ha,
                # cluster_rows = T,
                left_annotation = ha1,
                col=heat_colors2,
                color_space = "RGB",
                cluster_columns = FALSE,
                row_order=NULL,
                column_order= colnames(expp),
                show_column_names = T,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                gap = unit(1, "mm"),
                column_title_gp = gpar(fontsize = 7),
                width=unit(2.6, "cm"),
                height=unit(nrow(expp)*0.26, "cm"),
                show_heatmap_legend = T,
                heatmap_legend_param=list(labels_gp = gpar(fontsize = 7, fontface = "bold"),
                                          title_gp = gpar(fontsize = 7, fontface = "bold"))) 

anti
dev.off()

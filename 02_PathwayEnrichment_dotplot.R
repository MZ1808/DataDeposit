rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F)


if(!dir.exists('02_PathwayEnrichment_dotplot')){
  dir.create('02_PathwayEnrichment_dotplot')
}

index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SAE vs Ctrl.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

mytheme = theme(axis.text.x = element_text(hjust = 0.5,size = 8), 
                axis.text.y = element_text(size = 8), 
                axis.title.x = element_text(size = 8), 
                axis.title.y = element_text(size = 8), 
                # axis.line = element_line(size = 1),
                plot.margin = unit(c(1,1,1,1), "cm"),
                plot.title = element_text(hjust = 0.5,size = 8),
                legend.title = element_text(size = 8), 
                legend.text = element_text(size = 8), 
                legend.position = "right",
                legend.background = element_rect(fill = 'transparent'))

dt<-kegg
dt$GeneRatio<-c(as.numeric(dt$Count)/c(as.numeric(do.call(rbind,strsplit(dt$GeneRatio[1],'\\/'))[,2])))
dt$RichFactor <- dt$Count/as.numeric(sub("/\\d+", "", dt$BgRatio))
dt<-dt[order(dt$pvalue,decreasing = F),]
dt<-dt[1:20,]

dt2 <- dt[order(dt$RichFactor),]
dt2$Description <- factor(dt2$Description,levels = dt2$Description)            


p3 <- ggplot(data = dt2, aes(x = Description, y = RichFactor)) +
  geom_point(aes(color = pvalue, size = Count)) + 
  labs(x = NULL, y = "RichFactor", title = "KEGG Pathway Enrichment\nSAE vs Ctrl",color = "P value") +  
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3"),trans = "log10",
                        guide = guide_colorbar(reverse = TRUE,order = 1)) +
  scale_size_continuous(range = c(1, 5)) +                  
  coord_flip() + theme_bw() + mytheme +
  scale_x_discrete(labels = function(dat) str_wrap(dat,width = 30)) 
p3
ggsave(p3,filename = '02_PathwayEnrichment_dotplot/SAE vs Ctrl_dotplot.pdf',width = 5.3, height = 5.9)



index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SAE vs SP.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

dt<-kegg
dt$GeneRatio<-c(as.numeric(dt$Count)/c(as.numeric(do.call(rbind,strsplit(dt$GeneRatio[1],'\\/'))[,2])))
dt$RichFactor <- dt$Count/as.numeric(sub("/\\d+", "", dt$BgRatio))
dt<-dt[order(dt$pvalue,decreasing = F),]
dt<-dt[1:20,]

dt2 <- dt[order(dt$RichFactor),]
dt2$Description <- factor(dt2$Description,levels = dt2$Description)            


p3 <- ggplot(data = dt2, aes(x = Description, y = RichFactor)) +
  geom_point(aes(color = pvalue, size = Count)) +
  labs(x = NULL, y = "RichFactor", title = "KEGG Pathway Enrichment\nSAE vs SP",color = "P value") +  
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3"),trans = "log10",
                        guide = guide_colorbar(reverse = TRUE,order = 1)) +
  scale_size_continuous(range = c(1, 5)) +                  
  coord_flip() + theme_bw() + mytheme +
  scale_x_discrete(labels = function(dat) str_wrap(dat,width = 30)) 
p3
ggsave(p3,filename = '02_PathwayEnrichment_dotplot/SAE vs SP_dotplot.pdf',width = 5.3, height = 5.9)



index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SP vs Ctrl.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

dt<-kegg
dt$GeneRatio<-c(as.numeric(dt$Count)/c(as.numeric(do.call(rbind,strsplit(dt$GeneRatio[1],'\\/'))[,2])))
dt$RichFactor <- dt$Count/as.numeric(sub("/\\d+", "", dt$BgRatio))
dt<-dt[order(dt$pvalue,decreasing = F),]
dt<-dt[1:20,]

dt2 <- dt[order(dt$RichFactor),]
dt2$Description <- factor(dt2$Description,levels = dt2$Description)            


p3 <- ggplot(data = dt2, aes(x = Description, y = RichFactor)) +
  geom_point(aes(color = pvalue, size = Count)) + 
  labs(x = NULL, y = "RichFactor", title = "KEGG Pathway Enrichment\nSP vs Ctrl",color = "P value") +  
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3"),trans = "log10",
                        guide = guide_colorbar(reverse = TRUE,order = 1)) +
  scale_size_continuous(range = c(1, 5)) +                  
  coord_flip() + theme_bw() + mytheme +
  scale_x_discrete(labels = function(dat) str_wrap(dat,width = 30)) 
p3
ggsave(p3,filename = '02_PathwayEnrichment_dotplot/SP vs Ctrl_dotplot.pdf',width = 5.3, height = 5.9)

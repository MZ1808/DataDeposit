rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F)

if(!dir.exists('02_PathwayEnrichment_barplot')){
  dir.create('02_PathwayEnrichment_barplot')
}

index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SAE vs Ctrl.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

kegg$GeneRatio<-c(as.numeric(kegg$Count)/c(as.numeric(do.call(rbind,strsplit(kegg$GeneRatio[1],'\\/'))[,2])))
kegg$richFactor <- kegg$Count/as.numeric(sub("/\\d+", "", kegg$BgRatio))
kegg<-kegg[order(kegg$pvalue,decreasing = F),]
# kegg<-kegg[order(kegg$GeneRatio,decreasing = T),]
top20 <- data.frame(kegg$Description,kegg$richFactor ,kegg$pvalue)
top20<-top20[1:20,]
top20<-top20[order(top20$kegg.richFactor,decreasing = T),]
colnames(top20) <- c("Description","RichFactor","pval")

col_fun <- colorRampPalette(c("#ce342c","#fce0df"))(10) 
p <- ggplot(data=top20,aes(x=Description,y=RichFactor,fill=pval))

p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line.x = element_line(),
                 axis.line.y = element_line(),
                 legend.key.height = unit(2,'mm'),
                 legend.key.width = unit(3,'mm'),
                 legend.text = element_text(color="black",size=6),
                 legend.title = element_text(color="black",size=8),
                 axis.text=element_text(color="black",size=8),
                 axis.title.x = element_text(color="black",size=9)
)

p3 <- p2 + ylim(0,0.5) + scale_fill_gradientn(colours = col_fun)

p3 + scale_x_discrete(limits=rev(top20[,1])) +labs(x="",y="Rich factor",title="SAE vs Ctrl")
ggsave(filename = '02_PathwayEnrichment_barplot/SAE vs Ctrl_barplot.pdf',width = 6 , height = 3.53)


index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SAE vs SP.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

kegg$GeneRatio<-c(as.numeric(kegg$Count)/c(as.numeric(do.call(rbind,strsplit(kegg$GeneRatio[1],'\\/'))[,2])))
kegg$richFactor <- kegg$Count/as.numeric(sub("/\\d+", "", kegg$BgRatio))
kegg<-kegg[order(kegg$pvalue,decreasing = F),]
# kegg<-kegg[order(kegg$GeneRatio,decreasing = T),]

top20 <- data.frame(kegg$Description,kegg$richFactor ,kegg$pvalue)
top20<-top20[1:20,]
top20<-top20[order(top20$kegg.richFactor,decreasing = T),]
colnames(top20) <- c("Description","RichFactor","pval")

col_fun <- colorRampPalette(c("#ce342c","#fce0df"))(10)
p <- ggplot(data=top20,aes(x=Description,y=RichFactor,fill=pval))

p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line.x = element_line(),
                 axis.line.y = element_line(),
                 legend.key.height = unit(2,'mm'),
                 legend.key.width = unit(3,'mm'),
                 legend.text = element_text(color="black",size=6),
                 legend.title = element_text(color="black",size=8),
                 axis.text=element_text(color="black",size=8),
                 axis.title.x = element_text(color="black",size=9)
)

p3 <- p2 + ylim(0,0.8) + scale_fill_gradientn(colours = col_fun)

p3 + scale_x_discrete(limits=rev(top20[,1])) +labs(x="",y="Rich factor",title="SAE vs SP")
ggsave(filename = '02_PathwayEnrichment_barplot/SAE vs SP_barplot.pdf',width = 6.1 , height = 3.53)


index<-read_xlsx('02_PathwayEnrichment/enrichKEGG_SP vs Ctrl.xlsx',sheet = 1)
index<-as.data.frame(index)
# kegg<-data[data$ID%in%index$ID,]
kegg<-index
kegg<-kegg[kegg$pvalue<0.05,]
kegg$GeneRatio

kegg$GeneRatio<-c(as.numeric(kegg$Count)/c(as.numeric(do.call(rbind,strsplit(kegg$GeneRatio[1],'\\/'))[,2])))
kegg$richFactor <- kegg$Count/as.numeric(sub("/\\d+", "", kegg$BgRatio))
kegg<-kegg[order(kegg$pvalue,decreasing = F),]
# kegg<-kegg[order(kegg$GeneRatio,decreasing = T),]
top20 <- data.frame(kegg$Description,kegg$richFactor ,kegg$pvalue)
top20<-top20[1:20,]
top20<-top20[order(top20$kegg.richFactor,decreasing = T),]
colnames(top20) <- c("Description","RichFactor","pval")

col_fun <- colorRampPalette(c("#ce342c","#fce0df"))(10) 
p <- ggplot(data=top20,aes(x=Description,y=RichFactor,fill=pval))

#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line.x = element_line(),
                 axis.line.y = element_line(),
                 legend.key.height = unit(2,'mm'),
                 legend.key.width = unit(3,'mm'),
                 legend.text = element_text(color="black",size=6),
                 legend.title = element_text(color="black",size=8),
                 axis.text=element_text(color="black",size=8),
                 axis.title.x = element_text(color="black",size=9)
)

p3 <- p2 + ylim(0,0.4) + scale_fill_gradientn(colours = col_fun)

p3 + scale_x_discrete(limits=rev(top20[,1])) +labs(x="",y="Rich factor",title="SP vs Ctrl")
ggsave(filename = '02_PathwayEnrichment_barplot/SP vs Ctrl_barplot.pdf',width = 6.1 , height = 3.53)

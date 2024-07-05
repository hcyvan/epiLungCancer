source('R/base.R')
library(easyepi)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ComplexHeatmap)
library(patchwork)
library(circlize)
library(reshape2)

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
getRegionInLoyfer<-function(bed,split.category=TRUE){
  loyferT1000<-loadData(file.path(CONFIG$DataInter, 'dmc','loyfer','top1000_hg38.txt'),header = FALSE,ext='bed')
  gr1<-bed2gr(loyferT1000)
  gr2<-bed2gr(bed)
  overlaps<-findOverlaps(gr1,gr2)
  loyfer<-gr1[queryHits(overlaps)]
  lc<-gr2[subjectHits(overlaps)]
  num<-sort(table(loyfer$V4),decreasing = TRUE)
  category<-names(num)
  names(num)<-NULL
  data<-data.frame(
    category = category,
    value = as.vector(num)
  )
  if (split.category){
    data_split_category<-do.call(rbind,lapply(split(data, data$category), function(x){
      data.frame(
        category=strsplit(x$category,":")[[1]],
        value=x$value
      )
    }))
    data_split_category_merge<-do.call(rbind,lapply(split(data_split_category,data_split_category$category), function(x){
      data.frame(
        category=x$category[1],
        value=sum(x$value)
      )
    }))
    data<-data_split_category_merge
  }
  arrange(data,desc(value))
}

plotInterLoyfer<-function(data){
  data$fraction <- data$value / sum(data$value)
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$category, " ",   sprintf("%.1f%%", data$fraction * 100))
  
  ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
    scale_fill_brewer(palette=4) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void()+
    theme(legend.position = "none")
}

#'=============================================================================
#' Cell types and groups selected for deconvolution in this article
#'=============================================================================
typeSun<-c('Lungs','Brain','Heart','T-cells', 'B-cells','Neutrophils')
typeMoss<-c('LungCells','VascularEndothelialCells', 'LeftAtrium', 'CorticalNeurons','PancreaticBetaCells','Cd4tCells', 'Cd8tCells','BCells', 'NkCells','Monocytes', 'Neutrophils')
typeLoyfer<-c('Lung-Ep-Alveo', 'Lung-Ep-Bron', 'Endothel', 'Heart-Fibro', 'Head-Neck-Ep','Neuron', 'Oligodend', 'Pancreas-Alpha', 'Pancreas-Beta', 'Pancreas-Delta', 'Blood-T', 'Blood-B', 'Blood-NK', 'Blood-Mono+Macro', 'Blood-Granul')
typeGroup<-list(
  lung=list(
    sun=c('Lungs'),
    moss=c('LungCells'),
    loyfer=c('Lung-Ep-Alveo','Lung-Ep-Bron')
  ),
  neuroendocrine=list(
    sun=c('Brain'),
    moss=c('CorticalNeurons','PancreaticBetaCells'),
    loyfer=c('Neuron', 'Oligodend', 'Pancreas-Alpha', 'Pancreas-Beta', 'Pancreas-Delta')
  ),
  immune=list(
    sun=c('T-cells', 'B-cells','Neutrophils'),
    moss=c('Cd4tCells', 'Cd8tCells','BCells', 'NkCells','Monocytes', 'Neutrophils'),
    loyfer=c('Blood-T', 'Blood-B', 'Blood-NK', 'Blood-Mono+Macro', 'Blood-Granul')
  ),
  other=list(
    sun=c('Heart'),
    moss=c('VascularEndothelialCells','LeftAtrium'),
    loyfer=c('Endothel','Heart-Fibro', 'Head-Neck-Ep')
  )
)

typeGroupOrder<-c('lung','immune', 'neuroendocrine','other')
COLOR_MAP_TYPE<-c('#8ECFC9','#FFBE7A', '#FA7F6F','#82B0D2')
names(COLOR_MAP_TYPE)<-typeGroupOrder


getFraction<-function(method,mtype){
  data<-loadData(file.path(CONFIG$DataInter, 'deconv2',sprintf('deconv.%s.tsv',method)),header = TRUE,ext='bed')
  data$type<-factor(data$type, mtype)
  data<-arrange(data, type)
  m<-data[,-1]
  rownames(m)<-data$type
  m
}


DECOV<-list(
  sun=getFraction('sun',typeSun),
  moss=getFraction('moss',typeMoss),
  loyfer=getFraction('loyfer',typeLoyfer)
)

#'----------------------------------------------------------------------------------------------------------------------
#' Figure 3A. The overlaps of subtype-specific DMRs and cell type-specific Top1000 markers (Loyfer) 
#'----------------------------------------------------------------------------------------------------------------------
dmrHyperCTL<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.CTL.hyper.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLUAD<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LUAD.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLUSC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LUSC.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLCC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LCC.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoSCLC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.SCLC.hypo.dmr.bed'),header = TRUE,ext='bed')

dmrHyperCTL<-getRegionInLoyfer(dmrHyperCTL)
interHypoLUAD<-getRegionInLoyfer(dmrHypoLUAD)
interHypoLUSC<-getRegionInLoyfer(dmrHypoLUSC)
interHypoLCC<-getRegionInLoyfer(dmrHypoLCC)
interHypoSCLC<-getRegionInLoyfer(dmrHypoSCLC)

addFrac<-function(data){
  data$category<-factor(data$category, levels = data$category)
  data$fraction <- (data$value / sum(data$value))*100
  data$label <- sprintf("%.1f%%", data$fraction)
  data
}
interHypoLUAD<-addFrac(interHypoLUAD)
interHypoLUSC<-addFrac(interHypoLUSC)
interHypoLCC<-addFrac(interHypoLCC)
interHypoSCLC<-addFrac(interHypoSCLC)

interHypoLUAD<-interHypoLUAD[1:4,]
interHypoLUSC<-interHypoLUSC[1:4,]
interHypoLCC<-interHypoLCC[1:4,]
interHypoSCLC<-interHypoSCLC[1:4,]

interHypoLUAD$group<-'LUAD'
interHypoLUSC$group<-'LUSC'
interHypoLCC$group<-'LCC'
interHypoSCLC$group<-'SCLC'

interHypo<-do.call(rbind,list(
  interHypoLUAD,
  interHypoLUSC,
  interHypoLCC,
  interHypoSCLC
))
yl<-max(interHypo$fraction)*1.05
plotInterLoyferBarplot<-function(data,color="orange",show.xlab=TRUE){
  p<-ggplot(data, aes(x=category , y=fraction)) +
    geom_segment( aes(xend=category, yend=0)) +
    # geom_text(aes(y = fraction, x = category, label=label),vjust = -1, color = "black") +
    geom_point( size=4, color=color) +
    labs(x = NULL, y = NULL)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(0, yl))
  if (show.xlab){
    p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }else{
    p<-p+theme(axis.text.x = element_blank())
  }
  p
}

saveImage("loyfer.inter.barplot.hypo.pdf",width = 6,height = 1.3)
p1<-plotInterLoyferBarplot(interHypoLUAD,color = COLOR_MAP_GROUP[2])
p2<-plotInterLoyferBarplot(interHypoLUSC,color = COLOR_MAP_GROUP[3])
p3<-plotInterLoyferBarplot(interHypoLCC,color = COLOR_MAP_GROUP[4])
p4<-plotInterLoyferBarplot(interHypoSCLC,color = COLOR_MAP_GROUP[5])
grid.arrange(p1,p2,p3,p4, nrow = 1)
dev.off()

saveImage("loyfer.inter.barplot.hypo.without.xlab.pdf",width = 6,height = 1.3)
p1<-plotInterLoyferBarplot(interHypoLUAD,color = COLOR_MAP_GROUP[2],show.xlab = FALSE)
p2<-plotInterLoyferBarplot(interHypoLUSC,color = COLOR_MAP_GROUP[3],show.xlab = FALSE)
p3<-plotInterLoyferBarplot(interHypoLCC,color = COLOR_MAP_GROUP[4],show.xlab = FALSE)
p4<-plotInterLoyferBarplot(interHypoSCLC,color = COLOR_MAP_GROUP[5],show.xlab = FALSE)
grid.arrange(p1,p2,p3,p4, nrow = 1)
dev.off()

# saveImage("loyfer.inter.hypo.LUAD.pdf",width = 4,height = 3)
plotInterLoyfer(interHypoLUAD)
# dev.off()
#'----------------------------------------------------------------------------------------------------------------------
#' Figure S1A. Results of deconvolution of all samples using different algorithms
#'----------------------------------------------------------------------------------------------------------------------
plotHeatmap<-function(m, col_anno=FALSE){
  column_annotation<-NULL
  if (col_anno){
    column_annotation <-HeatmapAnnotation(
      df=data.frame(Stage=SAMPLE$sample2group(colnames(m))),
      col = list(Stage =COLOR_MAP_GROUP),
      show_annotation_name =FALSE,
      annotation_name_side='left'
    )
  }
  Heatmap(m,
          top_annotation = column_annotation,
          rect_gp = gpar(col = "black", lwd = 1),
          use_raster = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = FALSE
  )
}

p1<-plotHeatmap(DECOV$loyfer,TRUE)
p2<-plotHeatmap(DECOV$moss,TRUE)
p3<-plotHeatmap(DECOV$sun,TRUE)
saveImage("deconv.heatmap.pdf",width = 12,height = 6)
p1%v%p2%v%p3
dev.off()
#'----------------------------------------------------------------------------------------------------------------------
#' Figure S1B. Results of deconvolution of all samples using different algorithms
#'----------------------------------------------------------------------------------------------------------------------
getData<-function(data,method){
  dmerge<-sapply(names(COLOR_MAP_GROUP), function(x){
    dd<-SAMPLE$selectFromMatrixByGroup(data,x)
    rowMeans(dd)
  })
  m<-sapply(typeGroup, function(x){
    tgroup<-x[[method]]
    tmp<-data[match(tgroup,rownames(data)),]
    if (is.null(dim(tmp))){
      tmp
    }else{
      colSums(tmp)
    }
  })
  m
}

plotHeatmap<-function(m1,m2){
  p<-apply(m1, 2, function(col1) {
    apply(m2, 2, function(col2) {
      cor_result <- cor.test(col1, col2)
      cor_result$p.value
    })
  })
  
  corr<-apply(m1, 2, function(col1) {
    apply(m2, 2, function(col2) {
      cor_result <- cor.test(col1, col2)
      cor_result$estimate
    })
  })
  Heatmap(corr, name = "mat",
          col = colorRamp2(c(-1.2, 0, 1.2), c("green", "white", "red")),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_heatmap_legend =FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", corr[i, j]), x, y, gp = gpar(fontsize = 10))
          })
}
m1<-getData(DECOV$loyfer,'loyfer')
colnames(m1)<-paste0(colnames(m1),'1')
m2<-getData(DECOV$moss,'moss')
colnames(m2)<-paste0(colnames(m2),'2')
m3<-getData(DECOV$sun,'sun')
colnames(m3)<-paste0(colnames(m3),'3')
saveImage("deconv.heatmap.cor.1.pdf",width = 3.5,height = 3.5)
plotHeatmap(m1,m2)
dev.off()
saveImage("deconv.heatmap.cor.2.pdf",width = 3.5,height = 3.5)
plotHeatmap(m1,m3)
dev.off()
saveImage("deconv.heatmap.cor.3.pdf",width = 3.5,height = 3.5)
plotHeatmap(m2,m3)
dev.off()
#'----------------------------------------------------------------------------------------------------------------------
#' Figure 3B. Cellular components calculated by different deconvolution algorithms
#'----------------------------------------------------------------------------------------------------------------------
plotStackedBar<-function(data, method){
  dmerge<-sapply(names(COLOR_MAP_GROUP), function(x){
    dd<-SAMPLE$selectFromMatrixByGroup(data,x)
    rowMeans(dd)
  })
  m<-sapply(typeGroup, function(x){
    tgroup<-x[[method]]
    tmp<-dmerge[match(tgroup,rownames(dmerge)),]
    if (is.null(dim(tmp))){
      tmp
    }else{
      colSums(tmp)
    }
  })%>%t
  m<-data.frame(type=rownames(m),m)
  rownames(m)<-NULL
  data<-melt(m, variable.name = 'group')
  data$type <- factor(data$type, levels = typeGroupOrder)
  ggplot(data, aes(fill=type, y=value, x=group)) +
    scale_fill_manual(values = COLOR_MAP_TYPE)+
    geom_bar(position="fill", stat="identity")+
    labs(x = NULL, y = NULL)+
    theme_bw()+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}

saveImage("deconv.stacked.barplot.pdf",width = 6,height = 2)
p3<-plotStackedBar(DECOV$sun,'sun')
p2<-plotStackedBar(DECOV$moss,'moss')
p1<-plotStackedBar(DECOV$loyfer,'loyfer')
grid.arrange(p1,p2,p3, nrow = 1)
dev.off()

#'----------------------------------------------------------------------------------------------------------------------
#' Figure 3C. Boxplot of the cellular components of all samples
#'----------------------------------------------------------------------------------------------------------------------
plotCellTypePercentage<-function(data, ttype){
  frac<-data[match(ttype,rownames(data)),]
  frac<-unlist(frac)
  df<-data.frame(
    group=SAMPLE$sample2group(names(frac)),
    fraction=frac
  )
  ggplot(data=df,aes(x=group,y=fraction,fill=group))+
    scale_fill_manual(values=COLOR_MAP_GROUP) +
    geom_boxplot(outlier.shape = NA)+
    theme_bw()+
    labs(x="",y=sprintf('%s', ttype))+
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size=0),
          axis.title.y = element_text(size=14),
          axis.text.x = element_blank())+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}

plotBox<-function(data, types,nc,nr){
  combined_plot<-plotCellTypePercentage(data, types[1])
  for (ttype in types[2:length(types)]){
    p<-plotCellTypePercentage(data, ttype)
    combined_plot <- combined_plot + p
  }
  combined_plot <- combined_plot + plot_layout(ncol = nc, nrow = nr)
  print(combined_plot)
}

saveImage("deconv.boxplot.sun.pdf",width = 8,height = 6)
plotBox(DECOV$sun, typeSun, 4,4)
dev.off()
saveImage("deconv.boxplot.moss.pdf",width = 8,height = 6)
plotBox(DECOV$moss, typeMoss, 4,4)
dev.off()
saveImage("deconv.boxplot.loyfer.pdf",width = 8,height = 6)
plotBox(DECOV$loyfer, typeLoyfer, 4,4)
dev.off()

source('R/base.R')
library(easyepi)
library(ggplot2)
library(gridExtra)

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
plotInterLoyferBarplot<-function(data,color="orange"){
  data$category<-factor(data$category, levels = rev(data$category))
  data$fraction <- data$value / sum(data$value)
  data$label <- sprintf("%.1f%%", data$fraction * 100)
  ggplot(data, aes(x=category , y=value)) +
    geom_segment( aes(xend=category, yend=0)) +
    geom_text(aes(y = value, x = category, label=label),hjust = -0.2, color = "black") +
    geom_point( size=4, color=color) +
    coord_flip() +
    labs(x = "Cell Type", y = "Overlap Hypo Region") +
    theme_classic()+
    # theme_bw()+
    # theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(0, round(data$value[1]*1.1)))
}
#'----------------------------------------------------------------------------------------------------------------------
#' Figure 3A.
#'----------------------------------------------------------------------------------------------------------------------
dmrHyperCTL<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.CTL.hyper.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLUAD<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LUAD.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLUSC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LUSC.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoLCC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.LCC.hypo.dmr.bed'),header = TRUE,ext='bed')
dmrHypoSCLC<-loadData(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.SCLC.hypo.dmr.bed'),header = TRUE,ext='bed')

interHypoCTL<-getRegionInLoyfer(dmrHyperCTL)
interHypoLUAD<-getRegionInLoyfer(dmrHypoLUAD)
interHypoLUSC<-getRegionInLoyfer(dmrHypoLUSC)
interHypoLCC<-getRegionInLoyfer(dmrHypoLCC)
interHypoSCLC<-getRegionInLoyfer(dmrHypoSCLC)

interHypoLUAD<-interHypoLUAD[1:8,]
interHypoLUSC<-interHypoLUSC[1:8,]
interHypoLCC<-interHypoLCC[1:8,]
interHypoSCLC<-interHypoSCLC[1:8,]


saveImage("loyfer.inter.hypo.LUAD.pdf",width = 4,height = 3)
plotInterLoyfer(interHypoLUAD)
dev.off()

saveImage("loyfer.inter.barplot.hypo.pdf",width = 6,height = 6)
p1<-plotInterLoyferBarplot(interHypoLUAD,color = COLOR_MAP_GROUP[2])
p2<-plotInterLoyferBarplot(interHypoLUSC,color = COLOR_MAP_GROUP[3])
p3<-plotInterLoyferBarplot(interHypoLCC,color = COLOR_MAP_GROUP[4])
p4<-plotInterLoyferBarplot(interHypoSCLC,color = COLOR_MAP_GROUP[5])
grid.arrange(p1,p2,p3,p4, nrow = 2)
dev.off()





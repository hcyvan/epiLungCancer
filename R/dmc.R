source('R/base.R')
library(dplyr)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(easyepi)
#----------------------------------------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------------------------------------
getDmcFromMcomp<-function(mcomp.file, reverse=FALSE, has.percentile=FALSE){
  mcompDMC<-loadData(mcomp.file,header = TRUE,ext='bed')
  if (has.percentile){
    dmc<-dplyr::select(mcompDMC,chrom='#chrom',start,end, fc=`credibleDif_1-0`,p=`p_sim_1_v_0`,percentile,class)
  }else{
    dmc<-dplyr::select(mcompDMC,chrom='#chrom',start,end, fc=`credibleDif_1-0`,p=`p_sim_1_v_0`,class)
  }
  dmc$class<-ifelse(dmc$class=='strongHyper', 'hyper', ifelse(dmc$class=='strongHypo', 'hypo', 'nc'))
  if (reverse){
    dmc$fc<--dmc$fc
    dmc$class<-ifelse(dmc$class=='hyper', 'hypo', ifelse(dmc$class=='hypo', 'hyper', 'nc'))
  }
  return(dmc)
}
getDmcInfo<-function(dmc) {
  total<-sum(table(dmc$class))
  names(total)<-'total'
  return(c(table(dmc$class), total))
}
intersectBed<-function(x,y){
  inner_join(x, y, by = c('chrom'='chrom','start'='start', 'end'='end'))
}
multiInter<-function(bedList){
  sapply(names(bedList), function(x){
    sapply(names(bedList), function(y){
      intersectBed(bedList[[x]], bedList[[y]])%>%nrow
    })
  })
}
filterDmc<-function(dmcList,class=NULL,percentile=NULL, top=NULL){
  lapply(dmcList, function(x){
    out<-x
    if (!is.null(class)) {
      out<-out[out$class==class,]
    }
    if (!is.null(percentile)) {
      out<-out[out$percentile>=percentile,]
    }
    if (is.null(top)){
      return(out)
    }else{
      out<-arrange(out, desc(percentile),desc(abs(fc)))
      if (nrow(out)<=top) {
        return(out)
      }else {
        return(out[1:top,])
      }
    }
  })
}

multiInterAll<-function(dmcList){
  dmcListHyper<-filterDmc(dmcList, 'hyper')
  dmcListHypo<-filterDmc(dmcList, 'hypo')
  cat('Show all DMCs:\n')
  print(multiInter(dmcList))
  cat('Show Hyper DMCs:\n')
  print(multiInter(dmcListHyper))
  cat('Show Hypo DMCs:\n')
  print(multiInter(dmcListHypo))
}
removeMultiDMC<-function(bedList){
  cat('Covert bedList to longData ...\n')
  longData<-do.call(rbind,lapply(names(bedList), function(x){
    data<-bedList[[x]]
    data<-data.frame(feature=bed2Feature(data), group=x, value=1)
    return(data)
  }))
  cat('Covert longData to wideData ...\n')
  wideData<-dcast(longData,feature~group, value.var = 'value')
  cat('Remove duplicate interval occurrences ...\n')
  wideData[is.na(wideData)]<-0
  wideData2<-wideData[rowSums(wideData[,2:ncol(wideData)])==1,]
  cat('Covert wideData back to longData ...\n')
  longData2<-melt(wideData2,variable.name = "group")
  cat('Remove intervals that do not appear ...\n')
  longData2<-longData2[longData2$value==1,]
  cat('Covert longData to bed ...\n')
  final<-feature2Bed(longData2$feature)
  final2<-data.frame(final, group=longData2$group)
  final3<-split(final2, longData2$group)
  cat('Add information to the reserved interval ...\n')
  out<-lapply(names(bedList), function(x){
    tmp<-final3[[x]][1:3]
    return(left_join(tmp, bedList[[x]],by=c('chrom'='chrom','start'='start', 'end'='end')))
  })
  names(out)<-names(bedList)
  return(out)
}

plotPercentileVsKeepProp<-function(dmcWithPercentile,color='black',title=""){
  X<-seq(0, 1, by = 0.01)
  Y<-sapply(X, function(x){
    sum(dmcWithPercentile>=x)
  })
  
  df = data.frame(x=X*100, y=Y/Y[1])
  df<-df[df$x >=50,]
  ggplot(df, aes(x = x, y = y)) +
    geom_line(color = color) +
    labs(
      title = title,
      x = "Percentile",
      y = "Keep Proportion"
    )+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}
#'----------------------------------------------------------------------------------------------------------------------
#' Supplementary Material
#' Find the differentially methylated CpGs (DMCs)
#'----------------------------------------------------------------------------------------------------------------------
#'=============================================================================
#' The result of MOABS:mcomp
#'=============================================================================
one2rest<-list(
  RestvCTL=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.CTL.percentile.txt'),has.percentile = TRUE),
  RestvLUAD=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LUAD.percentile.txt'),has.percentile = TRUE),
  RestvLUSC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LUSC.percentile.txt'),has.percentile = TRUE),
  RestvLCC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LCC.percentile.txt'),has.percentile = TRUE),
  RestvSCLC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.SCLC.percentile.txt'),has.percentile = TRUE)
)
multiInterAll(one2rest)
#'=============================================================================
#' Remove DMCs appear in multi-group
#'=============================================================================
one2restRmMulti<-removeMultiDMC(one2rest)
multiInterAll(one2restRmMulti)
saveRDS(one2restRmMulti,file.path(CONFIG$DataInter, 'one2rest.0.rds'))
#'=============================================================================
#' Filter out DMCs by p-th percentile
#'=============================================================================
one2rest<-readRDS(file.path(CONFIG$DataInter, 'one2rest.0.rds'))
p1<-plotPercentileVsKeepProp(one2rest$RestvCTL$percentile, COLOR_MAP_GROUP[1],'CTL')
p2<-plotPercentileVsKeepProp(one2rest$RestvLUAD$percentile, COLOR_MAP_GROUP[2],'LUAD')
p3<-plotPercentileVsKeepProp(one2rest$RestvLUSC$percentile, COLOR_MAP_GROUP[3],'LUSC')
p4<-plotPercentileVsKeepProp(one2rest$RestvLCC$percentile, COLOR_MAP_GROUP[4],'LCC')
p5<-plotPercentileVsKeepProp(one2rest$RestvSCLC$percentile, COLOR_MAP_GROUP[5],'SCLC')
png(filename = "./notebook/img/dmcPercentileVsKeepProp.png", width = 8, height = 4, units = "in", res = 300)
grid.arrange(p1,p2,p3,p4,p5, nrow = 1)
dev.off()

one2rest<-readRDS(file.path(CONFIG$DataInter, 'one2rest.0.rds'))
one2rest<-filterDmc(one2rest, percentile = 0.8)
multiInterAll(one2rest)
saveRDS(one2rest,file.path(CONFIG$DataInter, 'one2rest.rds'))
#'=============================================================================
#' Generate a file containing DMCs for the extraction of CpGs methylation 
#' levels and the identification of differentially methylated regions (DMRs).
#' - methytools extract: extract of the CpGs methylation levels
#' - mcomppost dmc2dmr: merge DMCs into DMRs
#'=============================================================================
one2rest<-readRDS(file.path(CONFIG$DataInter, 'one2rest.rds'))
one2restDf<-do.call(rbind,lapply(names(one2rest), function(x){
  out<-one2rest[[x]]
  out$group<-x
  return(out)
}))
one2restDf$chrom<-factor(one2restDf$chrom,levels = HUMAN_CHROMSOME)
one2restDf<-arrange(one2restDf,chrom,start)
out<-one2restDf[,1:3]
out$class<-paste(one2restDf$class, one2restDf$group,sep = '.')
out$start <- format(out$start, scientific = FALSE)
out$end <- format(out$end, scientific = FALSE)
saveBed(out,file.path(CONFIG$DataInter, 'one2rest.dmc.bed'))
#'----------------------------------------------------------------------------------------------------------------------
#' Count of DMCs
#'----------------------------------------------------------------------------------------------------------------------
one2rest<-readRDS(file.path(CONFIG$DataInter, 'one2rest.rds'))
count<-sapply(one2rest, function(x){
  table(x$class)
})
count<-log10(count)
df <- data.frame(
  category = factor(c('CTL', 'LUAD','LUSC','LCC','SCLC','CTL', 'LUAD','LUSC','LCC','SCLC'),levels = FACTOR_LEVEL_GROUP),
  value = c(unlist(count[1,]),-unlist(count[2,])),
  group = rep(c("Hyper", "Hypo"), each = 5)
)
saveImage("dmc.count.pdf",width = 4,height = 3)
ggplot(df, aes(x = category, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(breaks = seq(-10, 10, 1)) +
  theme_minimal() +
  ylab("Log10(DMCs Number)") +
  scale_fill_manual(values = c("Hyper" = "#c82423", "Hypo" = "#2878b5"))+
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )
dev.off()

#'----------------------------------------------------------------------------------------------------------------------
#' DMC heatmap
#'----------------------------------------------------------------------------------------------------------------------
one2rest<-readRDS(file.path(CONFIG$DataInter, 'one2rest.rds'))


plotTopkHeatmap<-function(class='hypo', top=2000){
  hypo1000<-filterDmc(one2rest, class=class,top=top)
  multiInterAll(hypo1000)
  one2rest.dmc.beta<-loadData(file.path(CONFIG$DataInter, 'one2rest.dmc.beta.bed'),header = TRUE,ext='bed')
  features<-bed2feature(one2rest.dmc.beta)
  hypo1000M<-lapply(hypo1000, function(x){
    one2rest.dmc.beta[match(bed2feature(x),features),-c(1:3)]
  })
  m<-do.call(rbind,hypo1000M)
  m<-m[rowSums(m==-1)==0,]
  groups<-SAMPLE$sample2group(colnames(m))
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=groups),
    col = list(Stage =COLOR_MAP_GROUP),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  row_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=factor(sapply(strsplit(rownames(m),'\\.'),function(x){substr(x[1],6, nchar(x[1]))}),levels = FACTOR_LEVEL_GROUP)),
    col = list(Stage =COLOR_MAP_GROUP),
    show_annotation_name =FALSE,
    which = 'row'
  )
  
  m<-t(scale(t(m)))
  Heatmap(m,
          use_raster = FALSE,
          #col=colorRamp2(c(0, 0.5,1), c("#4574b6", "#fdfec2", "#d83127")),
          col=colorRamp2(c(-1, 0,1), c("#4574b6", "#fdfec2", "#d83127")),
          right_annotation  = row_annotation,
          top_annotation = column_annotation,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          # clustering_method_columns = "ward.D",
          show_row_names = FALSE,
          show_column_names = FALSE
  )
}
saveImage("dmc.heatmap.hypo.top2000.pdf",width = 5,height = 3.8)
plotTopkHeatmap('hypo',2000)
dev.off()
saveImage("dmc.heatmap.hyper.top1000.pdf",width = 5,height = 3.8)
plotTopkHeatmap('hyper',1000)
dev.off()
























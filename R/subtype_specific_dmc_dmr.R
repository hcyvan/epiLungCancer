source('R/base.R')
library(dplyr)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(easyepi)
library(rGREAT)
library(readxl)
# Function ---------------------------------------------------------------------
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
    data<-data.frame(feature=bed2feature(data), group=x, value=1)
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
  final<-feature2bed(longData2$feature)
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

annoDistQueryToSubject<-function(query, subject){
  midpoint<-function(gr){
    mid<-floor((start(gr) + end(gr)) / 2)
    start(gr)<-mid
    end(gr)<-mid
    gr
  }
  subject<-midpoint(subject)
  hits<-distanceToNearest(query, subject)
  query.hits<-query[queryHits(hits)]
  subject.hits<-subject[subjectHits(hits)]
  dist.sign<-ifelse((start(query.hits) > end(subject.hits))&(mcols(hits)$distance!=0),-1,1)
  query.hits$distanceToSubject<-mcols(hits)$distance * dist.sign
  query.hits
}
# Selection of subtype-specific DMCs -------------------------------------------
## The result of MOABS:mcomp and mcomppost percentile --------------------------
one2rest<-list(
  CTL=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.CTL.percentile.txt'),has.percentile = TRUE),
  LUAD=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LUAD.percentile.txt'),has.percentile = TRUE),
  LUSC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LUSC.percentile.txt'),has.percentile = TRUE),
  LCC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.LCC.percentile.txt'),has.percentile = TRUE),
  SCLC=getDmcFromMcomp(file.path(CONFIG$DataRaw,'moabs','dmc', 'dmc.Rest.vs.SCLC.percentile.txt'),has.percentile = TRUE)
)
multiInterAll(one2rest)
## Remove DMCs appear in multi-group and sex chromosomes -----------------------
one2restRmMulti<-removeMultiDMC(one2rest)
one2restRmMultiAndSex<-lapply(one2restRmMulti, function(x){
  filter(x, !chrom%in%c('chrX','chrY'))
})
multiInterAll(one2restRmMultiAndSex)
saveRDS(one2restRmMultiAndSex,file.path(CONFIG$DataInter, 'dmc','one2rest.rds'))
## Filter out DMCs by p-th percentile ------------------------------------------
one2rest<-readRDS(file.path(CONFIG$DataInter, 'dmc','one2rest.rds'))
p1<-plotPercentileVsKeepProp(one2rest$CTL$percentile, COLOR_MAP_GROUP[1],'CTL')
p2<-plotPercentileVsKeepProp(one2rest$LUAD$percentile, COLOR_MAP_GROUP[2],'LUAD')
p3<-plotPercentileVsKeepProp(one2rest$LUSC$percentile, COLOR_MAP_GROUP[3],'LUSC')
p4<-plotPercentileVsKeepProp(one2rest$LCC$percentile, COLOR_MAP_GROUP[4],'LCC')
p5<-plotPercentileVsKeepProp(one2rest$SCLC$percentile, COLOR_MAP_GROUP[5],'SCLC')
png(filename = "./notebook/img/dmcPercentileVsKeepProp.png", width = 8, height = 4, units = "in", res = 300)
grid.arrange(p1,p2,p3,p4,p5, nrow = 1)
dev.off()

one2rest<-readRDS(file.path(CONFIG$DataInter, 'dmc','one2rest.rds'))
multiInterAll(one2rest)
one2rest<-filterDmc(one2rest, percentile = 0.8)
multiInterAll(one2rest)
saveRDS(one2rest,file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.rds'))
one2rest<-filterDmc(one2rest, percentile = 0.85)
multiInterAll(one2rest)
saveRDS(one2rest,file.path(CONFIG$DataInter, 'dmc','p85','one2rest85.rds'))
one2rest<-filterDmc(one2rest, percentile = 0.90)
multiInterAll(one2rest)
saveRDS(one2rest,file.path(CONFIG$DataInter, 'dmc','p90','one2rest90.rds'))
# DMC to DMR -------------------------------------------------------------------
#'=============================================================================
#' Generate a file containing DMCs for the extraction of CpGs methylation 
#' levels and the identification of differentially methylated regions (DMRs).
#' - methytools extract: extract of the CpGs methylation levels
#' - mcomppost dmc2dmr: merge DMCs into DMRs
#'=============================================================================
outputDMC<-function(one2rest, outfile){
  one2restDf<-do.call(rbind,lapply(names(one2rest), function(x){
    out<-one2rest[[x]]
    out$group<-x
    return(out)
  }))
  one2restDf$chrom<-factor(one2restDf$chrom,levels = HUMAN_CHROMSOME)
  one2restDf<-arrange(one2restDf,chrom,start)
  out<-one2restDf[,1:3]
  out$class<-paste(one2restDf$class, one2restDf$group,sep = '.')
  out$start <- sapply(out$start, function(x){format(x, scientific = FALSE)})
  out$end <- sapply(out$end, function(x){format(x, scientific = FALSE)})
  out$percentile<-one2restDf$percentile
  saveBed(out,outfile)
}
one2rest<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.rds'))
outputDMC(one2rest,file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmc.bed'))
one2rest<-readRDS(file.path(CONFIG$DataInter, 'dmc','p85','one2rest85.rds'))
outputDMC(one2rest,file.path(CONFIG$DataInter, 'dmc','p85','one2rest85.dmc.bed'))
one2rest<-readRDS(file.path(CONFIG$DataInter, 'dmc','p90','one2rest90.rds'))
outputDMC(one2rest,file.path(CONFIG$DataInter, 'dmc','p90','one2rest90.dmc.bed'))
## workflow to get DMR ---------------------------------------------------------
#'------------------------------------------------------------------------------
#' workflow/methylation-level.smk
#'------------------------------------------------------------------------------
### mcomppost dmc2dmr ----------------------------------------------------------
#'------------------------------------------------------------------------------
#' mcomppost dmc2dmr -i one2rest80.dmc.bed -o one2rest80.dmr.py.bed
#'------------------------------------------------------------------------------
splitPyDmr<-function(dmrpy){
  dmr<-loadData(dmrpy,header = FALSE,ext='bed')
  dmr$group<-sapply(strsplit(dmr$V4,'\\.'),function(x){substr(x[2],6, nchar(x[2]))})
  dmr$class<-sapply(strsplit(dmr$V4,'\\.'),function(x){x[1]})
  dmr<-dplyr::select(dmr, chrom=V1, start=V2, end=V3, group, class, count=V5,length=V6)
  filename<-basename(dmrpy)
  outdir<-dirname(dmrpy)
  saveBed(dmr,file.path(outdir, gsub('.py.bed','.bed',filename)))
  dmrList<-list()
  for(g in FACTOR_LEVEL_GROUP){
    dmrSubList<-list()
    for(t in c('hypo','hyper')){
      dmr.tmp<-filter(dmr, class==t, group==g)
      file.name<-gsub('.dmr.py.bed',sprintf('.%s.%s.dmr.bed',g,t),filename)
      saveBed(dmr.tmp,file.path(outdir,file.name))
      dmrSubList[[t]]<-dmr.tmp
    }
    dmrList[[g]]<-dmrSubList
  }
  saveRDS(dmrList,file.path(outdir,gsub('.py.bed','.list.rds',filename)))
}
splitPyDmr(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.py.bed'))
splitPyDmr(file.path(CONFIG$DataInter, 'dmc','p85','one2rest85.dmr.py.bed'))
splitPyDmr(file.path(CONFIG$DataInter, 'dmc','p90','one2rest90.dmr.py.bed'))
### Top 500 subtype-specific DMR -----------------------------------------------
CTL.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.CTL.top500.hypo.dmr.bed"), header=FALSE, force.refresh = TRUE)
CTL.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.CTL.top500.hyper.dmr.bed"), header=FALSE, force.refresh = TRUE)
LUAD.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUAD.top500.hypo.dmr.bed"), header=FALSE, force.refresh = TRUE)
LUAD.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUAD.top500.hyper.dmr.bed"), header=FALSE, force.refresh = TRUE)
LUSC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUSC.top500.hypo.dmr.bed"), header=FALSE, force.refresh = TRUE)
LUSC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUSC.top500.hyper.dmr.bed"), header=FALSE, force.refresh = TRUE)
LCC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LCC.top500.hypo.dmr.bed"), header=FALSE, force.refresh = TRUE)
LCC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LCC.top500.hyper.dmr.bed"), header=FALSE, force.refresh = TRUE)
SCLC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.SCLC.top500.hypo.dmr.bed"), header=FALSE, force.refresh = TRUE)
SCLC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.SCLC.top500.hyper.dmr.bed"), header=FALSE, force.refresh = TRUE)

CTL.hypo$group='CTL'
LUAD.hypo$group='LUAD'
LUSC.hypo$group='LUSC'
LCC.hypo$group='LCC'
SCLC.hypo$group='SCLC'
CTL.hypo$group='CTL'
LUAD.hypo$group='LUAD'
LUSC.hypo$group='LUSC'
LCC.hypo$group='LCC'
SCLC.hypo$group='SCLC'

top500.dmr<-list(
  CTL=list(hypo=CTL.hypo,hyper=CTL.hyper),
  LUAD=list(hypo=LUAD.hypo,hyper=LUAD.hyper),
  LUSC=list(hypo=LUSC.hypo,hyper=LUSC.hyper),
  LCC=list(hypo=LCC.hypo,hyper=LCC.hyper),
  SCLC=list(hypo=SCLC.hypo,hyper=SCLC.hyper)
)
saveRDS(top500.dmr,file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.top500.list.rds'))
## DMR Density near TSS and CGI ------------------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$DataRaw,'genomicRegion.rds'))
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss

saveImage("dmr.density.genomicRegion.pdf",width = 5,height = 5.6)
par(mfrow = c(2, 2))
plotBedListDensity(dmrList,'hypo','tss', lim=5000, "Hypo",xlab='Distance to TSS')
plotBedListDensity(dmrList,'hypo','cgi', lim=5000,title="Hypo",xlab='Distance to CGI')
plotBedListDensity(dmrList,'hyper','tss', lim=5000, "Hyper",xlab='Distance to TSS')
plotBedListDensity(dmrList,'hyper','cgi', lim=5000, "Hyper",xlab='Distance to CGI')
dev.off()

# GREAT analysis ---------------------------------------------------------------
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.top500.list.rds'))
greatBP<-lapply(names(dmrList), function(x){
  lapply(names(dmrList[[x]]), function(y){
    dmr<-dmrList[[x]][[y]]
    print(paste(x, y))
    gr<-bed2gr(dmr)
    great(gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  })
})
names(greatBP)<-FACTOR_LEVEL_GROUP
for(g in FACTOR_LEVEL_GROUP){
  names(greatBP[[g]])<-c('hypo','hyper')
}
#saveRDS(greatBP,file.path(CONFIG$DataInter, 'dmc','p80','great','one2rest80.dmr.top500.list.greatBP.rds'))

greatBP<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','great','one2rest80.dmr.top500.list.greatBP.rds'))
res<-greatBP$CTL$hypo
getEnrichmentTable(res)
getRegionGeneAssociations(res)

greatBPtb<-do.call(rbind,lapply(names(greatBP), function(group){
  hypo<-greatBP[[group]][['hypo']]
  hypo<-getEnrichmentTable(hypo)
  hypo$class<-'hypo'
  hypo$group<-group
  hypo$key<-paste(group,'hypo',sep='.')
  hyper<-greatBP[[group]][['hyper']]
  hyper<-getEnrichmentTable(hyper)
  hyper$class<-'hyper'
  hyper$group<-group
  hyper$key<-paste(group,'hyper',sep='.')
  out<-rbind(hypo,hyper)
  out<-filter(out,p_adjust<=0.01,p_adjust_hyper<=0.01)
}))
write.csv(greatBPtb,file.path(CONFIG$DataInter,'dmc','p80','great','dmr.top500.greatBP.csv'),row.names = FALSE)
greatBPtbTop<-do.call(rbind,lapply(split(greatBPtb,greatBPtb$key), function(x){
  x[1:20,]
}))
write.csv(greatBPtbTop,file.path(CONFIG$DataInter,'dmc','p80','great','dmr.greatBP.top.csv'),row.names = FALSE)
# HOMER analysis ---------------------------------------------------------------
LUAD.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LUAD.top500.hypo.dmr'),q.value=0.05)
LUSC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LUSC.top500.hypo.dmr'),q.value=0.05)
LCC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LCC.top500.hypo.dmr'),q.value=0.05)
SCLC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.SCLC.top500.hypo.dmr'),q.value=0.05)
LUAD.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LUAD.top500.hyper.dmr'),q.value=0.05)
LUSC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LUSC.top500.hyper.dmr'),q.value=0.05)
LCC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.LCC.top500.hyper.dmr'),q.value=0.05)
SCLC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'/dmc/p80/homer/one2rest80.SCLC.top500.hyper.dmr'),q.value=0.05)

LUAD.hypo<-mutate(LUAD.hypo, group="LUAD",class="hypo")
LUAD.hyper<-mutate(LUAD.hyper, group="LUAD",class="hyper")
LUSC.hypo<-mutate(LUSC.hypo, group="LUSC",class="hypo")
LUSC.hyper<-mutate(LUSC.hyper, group="LUSC",class="hyper")
LCC.hypo<-mutate(LCC.hypo, group="LCC",class="hypo")
LCC.hyper<-mutate(LCC.hyper, group="LCC",class="hyper")
SCLC.hypo<-mutate(SCLC.hypo, group="SCLC",class="hypo")
SCLC.hyper<-mutate(SCLC.hyper, group="SCLC",class="hyper")

dmr.homer<-do.call(rbind,list(
  LUAD.hypo=LUAD.hypo,
  LUSC.hypo=LUSC.hypo,
  LCC.hypo=LCC.hypo,
  SCLC.hypo=SCLC.hypo,
  LUAD.hyper=LUAD.hyper,
  LUSC.hyper=LUSC.hyper,
  LCC.hyper=LCC.hyper,
  SCLC.hyper=SCLC.hyper
))
write.csv(dmr.homer,file.path(CONFIG$DataInter,'dmc','p80','homer','dmr.top500.homer80.csv'),row.names = FALSE)
# Figures ----------------------------------------------------------------------
## Figure 2B. Count of subtype-specific DMCs -----------------------------------
one2rest<-readRDS(file.path(CONFIG$DataInter,'dmc','p80','one2rest80.rds'))
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
## Figure 2C. subtype-specific DMCs Density near TSS and CGI --------
genomicRegion<-readRDS(file.path(CONFIG$DataRaw,'genomicRegion.rds'))
dmcList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.rds'))
dmcList<-lapply(dmcList, function(x){
  hypo<-filter(x,class=='hypo')
  hyper<-filter(x,class=='hyper')
  list(
    hypo=hypo,
    hyper=hyper
  )
})
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss

plot.dmc.density<-function(class,target.gr, lim,title="",xlab="Distance"){
  dist<-lapply(dmcList, function(x){
    bed<-x[[class]]
    distTo<-annoDistQueryToSubject(bed2gr(bed),target.gr)
    distTo
  })
  
  dens<-lapply(dist,function(x){
    d0<-x[abs(x$distanceToSubject)<lim,]
    if(length(d0$distanceToSubject)>0){
      density(d0$distanceToSubject)
    }
  })
  xlim<-range(sapply(dens, "[", "x"))
  ylim<-range(sapply(dens, "[", "y"))
  plot(NA, main = title,
       xlab = xlab,
       ylab = "DMRs Density",
       xlim=xlim,
       ylim=ylim)
  for(i in 1:length(COLOR_MAP_GROUP)) {
    lines(dens[[names(COLOR_MAP_GROUP)[i]]], lwd = 2,col=COLOR_MAP_GROUP[i])
  }
  legend("topleft", legend=names(COLOR_MAP_GROUP), fill=COLOR_MAP_GROUP, bty = "n")
}
saveImage("dmc.density.genomicRegion.pdf",width = 5,height = 5.6)
par(mfrow = c(2, 2))
plot.dmc.density('hypo',tss.gr, lim=5000, "Hypo",xlab='Distance to TSS')
plot.dmc.density('hypo',cgIslands.gr, lim=5000,title="Hypo",xlab='Distance to CGI')
plot.dmc.density('hyper',tss.gr, lim=5000, "Hyper",xlab='Distance to TSS')
plot.dmc.density('hyper',cgIslands.gr, lim=5000, "Hyper",xlab='Distance to CGI')
dev.off()
## Figure 2D. subtype-specific DMR heatmap -------------------------------------
CTL.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.CTL.top500.hypo.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
CTL.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.CTL.top500.hyper.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LUAD.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUAD.top500.hypo.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LUAD.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUAD.top500.hyper.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LUSC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUSC.top500.hypo.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LUSC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LUSC.top500.hyper.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LCC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LCC.top500.hypo.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
LCC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.LCC.top500.hyper.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
SCLC.hypo<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.SCLC.top500.hypo.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)
SCLC.hyper<-loadData(file.path(CONFIG$DataInter, "dmc/p80/one2rest80.SCLC.top500.hyper.dmr.sample.beta.bed"), header=TRUE, force.refresh = TRUE)

CTL.hypo$group='CTL'
CTL.hyper$group='CTL'
LUAD.hypo$group='LUAD'
LUAD.hyper$group='LUAD'
LUSC.hypo$group='LUSC'
LUSC.hyper$group='LUSC'
LCC.hypo$group='LCC'
LCC.hyper$group='LCC'
SCLC.hypo$group='SCLC'
SCLC.hyper$group='SCLC'

hypo.list<-list(LUAD.hypo, LUSC.hypo,LCC.hypo,SCLC.hypo)
hyper.list<-list(CTL.hyper, LUAD.hyper, LUSC.hyper,LCC.hyper,SCLC.hyper)
hypo<-do.call(rbind,hypo.list)
hyper<-do.call(rbind,hyper.list)


getHeatmap<-function(data){
  m<-data[,c(-1:-3)]
  m<-m[,-ncol(m)]
  groups<-SAMPLE$sample2group(colnames(m))
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=groups),
    col = list(Stage =COLOR_MAP_GROUP),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  row_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=factor(data$group,levels = FACTOR_LEVEL_GROUP)),
    col = list(Stage =COLOR_MAP_GROUP),
    show_annotation_name =FALSE,
    which = 'row'
  )
  
  Heatmap(m,
          use_raster = FALSE,
          col=colorRamp2(c(0, 0.5,1), c("#4574b6", "#fdfec2", "#d83127")),
          right_annotation  = row_annotation,
          top_annotation = column_annotation,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE
  )
}

h1<-getHeatmap(hypo)
h2<-getHeatmap(hyper)
saveImage("dmr.heatmap.top500.hypo.pdf",width = 5,height = 3.5)
# h1%v%h2
h1
dev.off()
## Figure 2E. Great analysis of subtype-specific DMRs --------------------------
data0<-read_xlsx(file.path(CONFIG$DataInter,'dmc','p80','great','dmr.greatBP.final.xlsx'))
data<-dplyr::select(data0, description=description, group=key,fc=fold_enrichment_hyper, p=p_adjust_hyper)
data<-as.data.frame(data)
data$nlogp<--log(data$p)
colorMap<-c('#2878b5','#c82423','#ffb15f','#fa66b3', '#925ee0','#2878b5','#c82423','#ffb15f','#fa66b3', '#925ee0')
names(colorMap)<-c('CTL.hypo','LUAD.hypo','LUSC.hypo','LCC.hypo','SCLC.hypo','CTL.hyper','LUAD.hyper','LUSC.hyper','LCC.hyper','SCLC.hyper')
data$group<-factor(data$group, levels = names(colorMap))
data$description<-factor(data$description, levels = rev(data$description))
saveImage("dmr.great.pdf",width = 7,height = 6)
ggplot(data,aes(x=group,y=description,size=nlogp,fill=group))+
  scale_fill_manual(values=colorMap)+
  geom_point(shape=21,color = "transparent")+
  theme(panel.background = element_blank(),
        panel.grid = element_line("gray"),
        panel.border = element_rect(colour = "black",fill=NA))
dev.off()

# Tables -----------------------------------------------------------------------
## Table S4. Subtype-specific DMR counts ---------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$DataRaw,'genomicRegion.rds'))
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
countGenomeDist<-function(bed){
  gr<-bed2gr(bed)
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr, gr1)!=0)
  })
}
dmrHypoCountInGenome<-sapply(names(dmrList), function(x){
  dmr<-dmrList[[x]][['hypo']]
  countGenomeDist(dmr)
})%>%t
dmrHyperCountInGenome<-sapply(names(dmrList), function(x){
  dmr<-dmrList[[x]][['hyper']]
  countGenomeDist(dmr)
})%>%t
dmrHypoCount<-sapply(names(dmrList), function(x){
  dmr<-dmrList[[x]][['hypo']]
  nrow(dmr)
})
dmrHyperCount<-sapply(names(dmrList), function(x){
  dmr<-dmrList[[x]][['hyper']]
  nrow(dmr)
})
dmrHypoCount<-data.frame(dmrHypoCount,dmrHypoCountInGenome)
dmrHyperCount<-data.frame(dmrHyperCount,dmrHyperCountInGenome)
write.csv(dmrHypoCount,file.path(CONFIG$DataResult,'table','dmr.hypo.count.csv'))
write.csv(dmrHyperCount,file.path(CONFIG$DataResult,'table','dmr.hyper.count.csv'))
## Table S5. Top 500 subtype-specific DMR  -------------------------------------
dmr.top500<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.top500.list.rds'))
LUAD.hypo<-dmr.top500$LUAD$hypo
LUSC.hypo<-dmr.top500$LUSC$hypo
LCC.hypo<-dmr.top500$LCC$hypo
SCLC.hypo<-dmr.top500$SCLC$hypo
LUAD.hypo$group='LUAD'
LUSC.hypo$group='LUSC'
LCC.hypo$group='LCC'
SCLC.hypo$group='SCLC'
hypo.list<-list(LUAD.hypo, LUSC.hypo,LCC.hypo,SCLC.hypo)
hypo<-do.call(rbind,hypo.list)
colnames(hypo)<-c('Chrom','Start','End','xx','CpGNumber','Length','Pt','Subtype')
write.csv(hypo,file.path(CONFIG$DataResult,'table','dmr.hypo.table5.hypo.top500.csv'))


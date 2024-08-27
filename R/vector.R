source('R/base.vector.R')


# Window size Attempts ---------------------------------------------------------
windowcount<-sapply(3:8, function(x){
  27852739-22*(x-1)
})
wins<-paste0('w',as.character(3:8))
wins<-factor(wins, levels = wins)
names(windowcount)<-wins
## Table S5 windows QC ---------------------------------------------------------
qcInfo<-sapply(wins, function(x){
  qc.file<-file.path(CONFIG$DataInter,'vector',x,'lung.mvc.qc')
  qc<-read.csv(qc.file, sep = '\t', header = TRUE, check.names = FALSE)
  colnames(qc)<-c('window','totalMVs', 'count','countNotEmpty','avgMVs','windowBpAvg')
  avg<-sum(as.numeric(strsplit(qc$total, split='\\|')[[1]]))/as.numeric(qc$countNotEmpty)
  cov<-qc$countNotEmpty/windowcount[match(x, names(windowcount))]
  ret<-c(qc$window, qc$windowBpAvg, avg, qc$countNotEmpty, cov)
  names(ret)<-c('CpGs','windowSize', 'MVs', 'windowNumber', 'Coverage')
  ret
})%>%t

qcInfo<-as.data.frame(qcInfo)
write.csv(qcInfo,file.path(CONFIG$DataResult,'table','mv.wins.qc.csv'))
## Figure 4C MVs number in windows ---------------------------------------------
avg<-lapply(wins, function(x){
  qc.file<-file.path(CONFIG$DataInter,'vector',x,'lung.mvc.qc')
  qc<-read.csv(qc.file, sep = '\t', header = TRUE, check.names = FALSE)
  colnames(qc)<-c('window','total', 'count','countNotEmpty','agv')
  sum(as.numeric(strsplit(qc$total, split='\\|')[[1]]))/as.numeric(qc$countNotEmpty)
})%>%unlist()
names(avg)<-3:8

lapply(wins, function(x){
  qc.file<-file.path(CONFIG$DataInter,'vector',x,'lung.mvc.qc')
  qc<-read.csv(qc.file, sep = '\t', header = TRUE, check.names = FALSE)
  colnames(qc)<-c('window','total', 'count','countNotEmpty','agv')
  qc$countNotEmpty
})%>%unlist()/windowcount

data <- data.frame(x=names(avg),y=avg)
saveImage("mv.average.MVs.pdf",width = 4,height = 3)
ggplot(data, aes(x=x, y=y)) +
  geom_segment(aes(x=x, xend=x, y=0, yend=y), color="grey") +
  geom_point(color="#b70e5e", size=4) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  xlab("CpGs count") +
  ylab("Average MVs Count")
dev.off()

qc.file<-file.path(CONFIG$DataInter,'vector/w6/lung.mvc.qc')
qc<-read.csv(qc.file, sep = '\t', header = TRUE, check.names = FALSE)
qc$agv
win4<-as.numeric(strsplit(qc$agv, split='\\|')[[1]])
pa<-plotMotifAnnoGgplot(6,'h')
motif<-Motif(6)
data<-data.frame(x=motif$get_motif_array(),y=win4)
p<-ggplot(data=data, aes(x=x, y=y)) +
  geom_bar(stat = "identity",fill='#b70e5e')+
  labs(x=NULL,y='MV counts')+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )

saveImage("mv.average.MVs.w4.pdf",width = 4.6,height = 2.6)
p/pa
dev.off()
## Figure S2 smvp=0.3 and ssp=0.9 -----------------------------------------------
plotFrac<-function(group="LUAD", window='w4'){
  vfile<-file.path(CONFIG$DataInter,'vector',window,paste0(group,'_0.80_0.1.mvc'))
  vector<-loadData(vfile,ext = 'bed')
  vector<-vector[,c(1,3,4,2,5:ncol(vector))]
  vector<-(filter(vector, V11>=0.9, V12>0.3))
  data<-data.frame(
    x=vector$V11,
    y=vector$V12
  )
  p_scatter <- ggplot(data, aes(x = x, y = y)) +
    geom_point(size=0.3,alpha = 0.6) +
    labs(x = NULL, y = NULL) +
    ggtitle(nrow(data))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(
      plot.title = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12)
    )
  p<-ggMarginal(p_scatter, type="histogram",fill='white',bins = 100,size = 8)
  p
}

p1.3<-plotFrac('LUAD', 'w3')
p1.4<-plotFrac('LUAD', 'w4')
p1.5<-plotFrac('LUAD', 'w5')
p1.6<-plotFrac('LUAD', 'w6')
p1.7<-plotFrac('LUAD', 'w7')
p1.8<-plotFrac('LUAD', 'w8')

p2.3<-plotFrac('LUSC', 'w3')
p2.4<-plotFrac('LUSC', 'w4')
p2.5<-plotFrac('LUSC', 'w5')
p2.6<-plotFrac('LUSC', 'w6')
p2.7<-plotFrac('LUSC', 'w7')
p2.8<-plotFrac('LUSC', 'w8')

p3.3<-plotFrac('LCC', 'w3')
p3.4<-plotFrac('LCC', 'w4')
p3.5<-plotFrac('LCC', 'w5')
p3.6<-plotFrac('LCC', 'w6')
p3.7<-plotFrac('LCC', 'w7')
p3.8<-plotFrac('LCC', 'w8')

p4.3<-plotFrac('SCLC', 'w3')
p4.4<-plotFrac('SCLC', 'w4')
p4.5<-plotFrac('SCLC', 'w5')
p4.6<-plotFrac('SCLC', 'w6')
p4.7<-plotFrac('SCLC', 'w7')
p4.8<-plotFrac('SCLC', 'w8')

pp<-grid.arrange(
  p1.3,p1.4,p1.5,p1.6,p1.7,p1.8,
  p2.3,p2.4,p2.5,p2.6,p2.7,p2.8,
  p3.3,p3.4,p3.5,p3.6,p3.7,p3.8,
  p4.3,p4.4,p4.5,p4.6,p4.7,p4.8,
  ncol=6)
ggsave(file.path(CONFIG$DataResult, 'mv.smvp.ssp.png'), plot = pp, width = 8, height = 7, dpi = 300)
# saveImage("mv.smvp.ssp.pdf",width = 8,height = 7)
# pp
# dev.off()


vfile<-file.path(CONFIG$DataInter,'vector/w4',paste0('LUAD','_0.80_0.1.smvc'))
vector<-loadData(vfile,ext = 'bed')
vector<-vector[,c(1,3,4,2,5:ncol(vector))]
vector<-(filter(vector, V11>=0.9, V12>0.4))
dim(vector)
# SMVC -------------------------------------------------------------------
groups<-names(COLOR_MAP_GROUP)[-1]
files<-lapply(groups, function(x){
  list(
    a=file.path(CONFIG$DataInter,'vector/w4',paste0(x,'_1.0_0.4.smvc')),
    b=file.path(CONFIG$DataInter,'vector/w4',paste0(x,'.smvc'))
  )
})
# names(files)<-groups
# files[[3]]$a<-file.path(CONFIG$DataInter,'vector/w4','LCC_1.0_0.3.smvc')
smvcs<-lapply(files, function(x){
  SMVC(x$a,x$b)
})
saveRDS(smvcs,file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
## Save SMVC -------------------------------------------------------------------
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
._<-lapply(FACTOR_LEVEL_GROUP[-1], function(x){
  dmr<-dmrList[[x]]
  smvc<-smvcs[[x]]
  smvc$freeze(file.path(CONFIG$DataInter,'vector/w4/bed'), dmr)
})
## Intersect with DMR ----------------------------------------------------------
do.call(rbind,lapply(FACTOR_LEVEL_GROUP[-1], function(x){
  dmr<-dmrList[[x]]
  smvc<-smvcs[[x]]
  ov1<-smvc$overlap(smvc$hypo, dmr$hypo)
  ov2<-smvc$overlap(smvc$hyper, dmr$hyper)
  data<-rbind(c(ov1$gr0stat, ov1$gr1stat),c(ov2$gr0stat, ov2$gr1stat))
  colnames(data)<-c('mvc', 'mvcUniq', 'mvcInter','dmr','dmrUniq','dmrInter')
  data<-data.frame(group=x,data,class=c('hypo','hyper'))
  data
}))

## Table S7. Annotation -----------------------------------------------------------------
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
anno<-lapply(smvcs, function(x){
  annoPeak(bed2gr(x$table))
})

smvc.anno<-lapply(names(smvcs), function(x){
  smvc<-smvcs[[x]]
  f1<-bed2feature(smvc$table)
  an<-anno[[x]]
  an<-gr2bed(an)
  f2<-bed2feature(an)
  an<-an[match(f1,f2),]
  annotation<-an$annotation
  distanceToTSS<-an$distanceToTSS
  symbol<-an$symbol
  ensembl<-an$ensembl
  motif1<-read.csv(file.path(CONFIG$DataInter,'vector/w4/homer',sprintf('%s_1.0_0.4.hypo.gr0.txt',x)),sep = '\t')
  motif2<-read.csv(file.path(CONFIG$DataInter,'vector/w4/homer',sprintf('%s_1.0_0.4.hyper.gr0.txt',x)),sep = '\t')
  motifs<-rbind(motif1, motif2)
  tfmap<-sapply(split(motifs,motifs$PositionID), function(x){
    aa<-sapply(strsplit(x$Motif.Name,'\\('),function(x){x[1]})
    paste(unique(aa),collapse = ',')
  })%>%unlist
  tfs<-tfmap[match(f1,names(tfmap))]
  if (is.null(tfs)){
    tfs<-NA
  }
  out<-data.frame(
    dplyr::select(smvc$table, chrom, start,end,cpg, mvs, smvp,ssp, i0,i1, i0v0,i1v0,class),
    group=x,
    symbol=symbol,
    ensembl=ensembl,
    distanceToTSS=distanceToTSS,
    annotation=annotation,
    TFs=tfs
  )
  #data.frame(smvc$table, annotation, distanceToTSS, symbol,ensembl,tfs)
})

smvc.anno<-do.call(rbind,smvc.anno)
saveBed(smvc.anno,file.path(CONFIG$DataInter,'vector/w4/anno/smvc.anno.bed'))
## Visualization ---------------------------------------------------------------
### LUAD -----------------------------------------------------------------------
smvc<-SMVC(file.path(CONFIG$DataInter,'vector/w4/LUAD_1.0_0.4.smvc'))
beta<-Beta(file.path(CONFIG$DataInter,'vector/w4/beta/LUAD_1.0_0.4.sample.beta.bed'))
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/mvm/LUAD.sample.mvm'))

arrange(smvc$hypo, desc(ssp))%>%head
arrange(smvc$hyper, desc(ssp))%>%head

cpg<-2188467
bvalue<-beta$getBeta(smvc$cpg2genome(cpg))
data<-mvm$getByCpG(cpg)

#saveImage("mv.window.pdf",width = 4.6,height = 6)
plotWindow(data,4,bvalue)
#dev.off()

cpg<-9589060
bvalue<-beta$getBeta(smvc$cpg2genome(cpg))
data<-mvm$getByCpG(cpg)

#saveImage("mv.window.pdf",width = 4.6,height = 6)
plotWindow(data,4,bvalue)
#dev.off()
### SCLC -----------------------------------------------------------------------
#### Figure 5B,C ---------------------------------------------------------------
smvc<-SMVC(file.path(CONFIG$DataInter,'vector/w4/SCLC_1.0_0.4.smvc'))
beta<-Beta(file.path(CONFIG$DataInter,'vector/w4/beta/SCLC_1.0_0.4.sample.beta.bed'))
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/mvm/SCLC.sample.mvm'))
arrange(smvc$hypo, desc(ssp))%>%head #chr1:4505825
arrange(smvc$hyper, desc(ssp))%>%head #chr10:16373819
cpg<-4505825
cpg<-20775751
bvalue<-beta$getBeta(smvc$cpg2genome(cpg))
data<-mvm$getByCpG(cpg)
saveImage("smvc.SCLC.hypo.20775751.pdf",width = 5,height = 6)
plotWindow(data,4,bvalue)
dev.off()

cpg<-16373819
bvalue<-beta$getBeta(smvc$cpg2genome(cpg))
data<-mvm$getByCpG(cpg)
saveImage("smvc.SCLC.hyper.16373819.pdf",width = 5,height = 6)
plotWindow(data,4,bvalue)
dev.off()

#### Figure 5A smvc heatmap ----------------------------------------------------
beta<-Beta(file.path(CONFIG$DataInter,'vector/w4/beta/SCLC_1.0_0.4.group.beta.bed'))
mark<-c(4505825,20775751, 16373819)
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/mvm/SCLC.group.mvm'))
p0<-plotSmvcHeatmap(mvm$getBySample('CTL',4),4,'CTL',mark,beta$table$CTL)
p1<-plotSmvcHeatmap(mvm$getBySample('LUAD',4),4,'LUAD',mark,beta$table$LUAD)
p2<-plotSmvcHeatmap(mvm$getBySample('LUSC',4),4,'LUSC',mark,beta$table$LUSQ)
p3<-plotSmvcHeatmap(mvm$getBySample('LCC',4),4,'LCC',mark,beta$table$LCC)
p4<-plotSmvcHeatmap(mvm$getBySample('SCLC',4),4,'SCLC',mark,beta$table$SCLC)
saveImage("smvc.SCLC.pdf",width =12,height = 6)
p0%v%p1%v%p2%v%p3%v%p4
dev.off()

#### Figure 5D smvc region ----------------------------------------------------
#mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.4505825.chr1_4505805_4505845.group.mvm'))
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.16373819.chr10_16373799_16373839.group.mvm')) # OK
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.14745014.chr9_14744964_14745064.group.mvm')) # NeuroD1
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.14745014.chr9_147449684_14745044.group.mvm')) # NeuroD1
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.20775751.chr14_20775721_20775781.group.mvm')) # NeuroD1
#mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.4505825.chr1_4505775_4505875.group.mvm')) # ok

pp0<-plotSmvcRegion(mvm$getBySample('CTL'))
pp1<-plotSmvcRegion(mvm$getBySample('LUAD'))
pp2<-plotSmvcRegion(mvm$getBySample('LUSC'))
pp3<-plotSmvcRegion(mvm$getBySample('LCC'))
pp4<-plotSmvcRegion(mvm$getBySample('SCLC'))
saveImage("smvc.SCLC.hypo.20775751.chr14_20775721_20775781.pdf",width = 9,height = 7)
pp0/pp1/pp2/pp3/pp4
dev.off()

# HOMER analysis ---------------------------------------------------------------
fc<-1.3
q<-0.01

LUAD.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUAD_1.0_0.4.hypo.gr0.reduce'),q, fc)
LUAD.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUAD_1.0_0.4.hyper.gr0.reduce'),q, fc)
LUSC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUSC_1.0_0.4.hypo.gr0.reduce'),q, fc)
LUSC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUSC_1.0_0.4.hyper.gr0.reduce'),q, fc)
LCC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LCC_1.0_0.4.hypo.gr0.reduce'),q, fc)
LCC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LCC_1.0_0.4.hyper.gr0.reduce'),q, fc)
SCLC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/SCLC_1.0_0.4.hypo.gr0.reduce'),q, fc)
SCLC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/SCLC_1.0_0.4.hyper.gr0.reduce'),q, fc)

fmg<-list(
  LUAD.hypo=LUAD.hypo,
  LUSC.hypo=LUSC.hypo,
  LCC.hypo=LCC.hypo,
  SCLC.hypo=SCLC.hypo,
  LUAD.hyper=LUAD.hyper,
  LUSC.hyper=LUSC.hyper,
  LCC.hyper=LCC.hyper,
  SCLC.hyper=SCLC.hyper
)

output<-lapply(names(fmg), function(n){
  x<-fmg[[n]]
  if(nrow(x)>30){
    out<-x[1:30,]
  }else{
    out<-x
  }
  if (nrow(out)>0){
    out$class<-n
  }
  out
})
output<-do.call(rbind, output)
#write.csv(output,file.path(CONFIG$DataInter,'vector/w4/homer','smvc.homer.csv'),row.names = FALSE)
write.csv(output,file.path(CONFIG$DataInter,'vector/w4/homer','smvc.homer.reduce.csv'),row.names = FALSE)
## Figure 5E TFs insercet with DMR ----------------------------------------------------------

homerDMR<-read.csv(file.path(CONFIG$DataInter,'dmc','p80','homer','dmr.homer80.top.csv'))
homerSMVC<-read.csv(file.path(CONFIG$DataInter,'vector/w4/homer','smvc.homer.reduce.csv'))
homerDMR<-split(homerDMR, homerDMR$class)
homerSMVC<-split(homerSMVC, homerSMVC$class)

r1<-homerSMVC$LUSC.hypo
r2<-homerDMR$LUSC.hypo

statInter<-function(r1,r2){
  tf1<-unique(r1$tf)
  tf2<-unique(r2$tf)
  
  tf1tf2<-intersect(tf1,tf2)
  tf1u<-setdiff(tf1,tf1tf2)
  tf2u<-setdiff(tf2,tf1tf2)
  
  cat(sprintf("%d\t%d\t%d\n", length(tf1u),length(tf1tf2), length(tf2u)))
  print(tf1u)
  print(tf1tf2)
  print(tf2u)
}

statInter(homerDMR$LUSC.hypo,homerSMVC$LUSC.hypo)
statInter(homerDMR$SCLC.hypo,homerSMVC$SCLC.hypo)
statInter(homerDMR$SCLC.hyper,homerSMVC$SCLC.hyper)
## Motif binding region --------------------------------------------------------
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
motif<-read.csv(file.path(CONFIG$DataInter,'vector/w4/homer/SCLC_1.0_0.4.hypo.gr0.reduce.txt'),sep='\t',check.names = FALSE)
motif$tf<-sapply(strsplit(motif$`Motif Name`,'\\('),function(x){x[1]})
motif<-filter(motif, tf=='NeuroD1')

bed<-smvcs$SCLC$reduceBed$hypo
bed$feature<-bed2feature(bed)

bed[match(unique(motif$PositionID),bed$feature),]
# 449,450,451,452  chr9:124330599-124330885
smvcs$SCLC$hypo[c(449,450,451,452),]

# SMVC Density near TSS and CGI -------------------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$DataRaw,'genomicRegion.rds'))
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
#saveImage("smvc.density.genomicRegion.pdf",width = 5,height = 5.6)
par(mfrow = c(2, 2))
plotBedListDensity(smvcs,'hypo','tss', lim=5000, "Hypo",xlab='Distance to TSS',COLOR_MAP_GROUP[-1])
plotBedListDensity(smvcs,'hypo','cgi', lim=5000, "Hypo",xlab='Distance to TSS',COLOR_MAP_GROUP[-1])
plotBedListDensity(smvcs,'hyper','tss', lim=5000, "Hyper",xlab='Distance to TSS',COLOR_MAP_GROUP[-1])
plotBedListDensity(smvcs,'hyper','cgi', lim=5000, "Hyper",xlab='Distance to TSS',COLOR_MAP_GROUP[-1])
dev.off()
# GREAT analysis ---------------------------------------------------------------
#dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.rds'))
## smvc -------------------------
greatBP<-lapply(names(smvcs), function(x){
  lapply(c('hypo','hyper'), function(y){
    dmr<-smvcs[[x]][[y]]
    print(paste(x, y))
    gr<-bed2gr(dmr)
    great(gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  })
})
names(greatBP)<-FACTOR_LEVEL_GROUP[-1]
for(g in FACTOR_LEVEL_GROUP[-1]){
  names(greatBP[[g]])<-c('hypo','hyper')
}
saveRDS(greatBP,file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.rds'))
## reduced smvc ----------------------------------------------------------------
greatBP.reduce<-lapply(names(smvcs), function(x){
  lapply(c('hypo','hyper'), function(y){
    gr<-smvcs[[x]]$reduceGr[[y]]
    print(paste(x, y))
    great(gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  })
})
names(greatBP.reduce)<-FACTOR_LEVEL_GROUP[-1]
for(g in FACTOR_LEVEL_GROUP[-1]){
  names(greatBP.reduce[[g]])<-c('hypo','hyper')
}
saveRDS(greatBP.reduce,file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.reduce.rds'))
## tmp --------------------------------------

saveGreatBP<-function(rds, csv,csv.top){
  greatBP<-readRDS(rds)
  greatBPtb<-do.call(rbind,lapply(names(greatBP), function(group){
    hypo<-greatBP[[group]][['hypo']]
    hypo<-getEnrichmentTable(hypo)
    hypo$class<-'hypo'
    hypo$group<-group
    hypo$key<-paste(group,'hypo',sep='.')
    hyper<-greatBP[[group]][['hyper']]
    hyper<-getEnrichmentTable(hyper)
    if(nrow(hyper)!=0){
      hyper$class<-'hyper'
      hyper$group<-group
      hyper$key<-paste(group,'hyper',sep='.')
      out<-rbind(hypo,hyper)
    }else{
      out<-hypo
    }
    out<-filter(out,p_adjust<=0.01,p_adjust_hyper<=0.01)
  }))
  write.csv(greatBPtb,csv,row.names = FALSE)
  greatBPtbTop<-do.call(rbind,lapply(split(greatBPtb,greatBPtb$key), function(x){
    x[1:20,]
  }))
  write.csv(greatBPtbTop,csv.top,row.names = FALSE)
}

saveGreatBP(file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.rds'),
            file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.csv'),
            file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.top.great.BP.csv'))

saveGreatBP(file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.reduce.rds'),
            file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.great.BP.reduce.csv'),
            file.path(CONFIG$DataInter, 'vector/w4/great/smvcs.top.great.BP.reduce.csv'))







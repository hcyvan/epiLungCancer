source('R/base.R')
library(easyepi)
library(openxlsx)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(rGREAT)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

## Function --------------------------------------------------------------------
sheets<-c(
  'finn-b-C3_SCLC_EXALLC',
  'finn-b-C3_LUNG_NEUROENDO_EXALLC',
  'finn-b-C3_NSCLC_SQUAM_EXALLC',
  'ieu-a-967',
  'ieu-a-989',
  'finn-b-C3_NSCLC_ADENO_EXALLC',
  'finn-b-C3_LUNG_NONSMALL_EXALLC',
  'ieu-a-966',
  'ieu-a-985',
  'ieu-a-986',
  'ieu-a-987',
  'ieu-b-4954',
  'ukb-a-54',
  'bbj-a-133'
)
# SMR results ------------------------------------------------------------------
anno<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
mqtls<-lapply(sheets, function(x){
  print(x)
  mqtl<-read.xlsx(file.path(CONFIG$DataInter,'smr','mqtl_smr_results.xlsx'),sheet = x)
  am<-anno[match(mqtl$probeID, rownames(anno)),]
  mqtl2<-data.frame(chrom=am$CHR_hg38,start=as.numeric(am$Start_hg38)-1, end=as.numeric(am$End_hg38),mqtl)
  mqtl2<-mqtl2[!is.na(mqtl2$chrom),]
  mqtl2<-filter(mqtl2,p_SMR < 0.05, p_HEIDI>0.05)
  mqtl2
})
names(mqtls)<-sheets
saveRDS(mqtls,file.path(CONFIG$Result,'mqtls.rds'))
## SRM statistics --------------------------------------------------------------
sapply(mqtls, function(x){
  nrow(x)
})
cat(sprintf('%d SMR results in %d datasets', nrow(mqtlTable), length(mqtls)))
cat(sprintf('%d CpGs in %d datasets', nrow(mqtlAll), length(mqtls))) #22303 


m<-sapply(mqtls, function(x){
  sapply(mqtls, function(y){
    length(intersect(x$probeID,y$probeID))
  })
})
saveImage("smr.numer.pdf",width = 6,height = 5)
Heatmap(m,
        col=colorRamp2(c(0, max(m)),c("#fffeee", "#c82423"))
        )
dev.off()

## Function Analysis -----------------------------------------------------------
### Merge adjacent cpgs --------------------------------------------------------
mqtls<-readRDS(mqtls,file.path(CONFIG$Result,'mqtls.rds'))
mqtlsMerged<-lapply(names(mqtls), function(x){
  mqtl<-mqtls[[x]]
  gr<-bed2gr(mqtl)
  gr<-reduce(gr,  min.gapwidth = 400)
  bed<-gr2bed(gr)
  saveBed(bed,file.path(CONFIG$DataInter, 'smr',sprintf('smr.%s.bed',x)))
  bed
})
names(mqtlsMerged)<-names(mqtls)
### rGREAT ---------------------------------------------------------------------
mqtlGreat400<-lapply(mqtlsMerged, function(x){
  gr<-bed2gr(x)
  list(
    bp=great(gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
    mf=great(gr, "GO:MF", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
    cc=great(gr, "GO:CC", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  )
})
saveRDS(mqtlGreat400,file.path(CONFIG$DataInter, 'smr','mqtlsGreat.g400.rds'))
mqtlGreat400<-saveRDS(file.path(CONFIG$DataInter, 'smr','mqtlsGreat.g400.rds'))

mqtlGreat400Filtered<-lapply(mqtlGreat400, function(x){
  lapply(c('bp','mf','cc'), function(y){
    aa<-getEnrichmentTable(x[[y]])
    aa<-filter(aa,p_adjust<=0.01,p_adjust_hyper<=0.01)
  })[[1]]
})
saveRDS(mqtlGreat400Filtered,file.path(CONFIG$Result,'mqtls.great.rds'))
mqtlGreat400Filtered<-readRDS(file.path(CONFIG$Result,'mqtls.great.rds'))

mqtlsGreatTable<-do.call(rbind,lapply(names(mqtlGreat400Filtered), function(x){
  g<-mqtlGreat400Filtered[[x]]
  data.frame(id=g$id,description=g$description, dataset=x, num=1)
}))

table(mqtlsGreatTable$description)%>%sort
table(mqtlsGreatTable$dataset)%>%sort

a<-dcast(mqtlsGreatTable,description~dataset,value.var = "num",fill=0)
m<-as.matrix(a[,-1])
rownames(m)<-a$description
m<-m[,match(sheets, colnames(m))]
a<-dcast(mqtlsGreatTable,id~dataset,value.var = "num",fill=0)
m2<-as.matrix(a[,-1])
rownames(m2)<-a$id
m2<-m2[,match(sheets, colnames(m2))]

saveImage("smr.great.3.pdf",width = 6,height = 30)
Heatmap(m,
        col=colorRamp2(c(0, 1),c("#fffeee", "#c82423")),
        cluster_columns = FALSE
)
dev.off()
saveImage("smr.great.2.pdf",width = 6,height = 30)
Heatmap(m2,
        col=colorRamp2(c(0, 1),c("#fffeee", "#c82423")),
        cluster_columns = FALSE
)
dev.off()
saveImage("smr.great.pdf",width = 5,height = 8)
Heatmap(m2,
        col=colorRamp2(c(0, 1),c("#fffeee", "#c82423")),
        cluster_columns = FALSE
)
dev.off()
### Homer ----------------------------------------------------------------------
fc<-0
p<-0.05
mqtlsHomer<-do.call(rbind,lapply(names(mqtlsMerged), function(x){
  homerResult<-file.path(CONFIG$DataInter, 'smr/homer',sprintf('smr.%s.bed',x))
  res<-getfindMotifsGenomeResults(homerResult,p.value=p, foldchange = fc)
  tf<-lapply(strsplit(res$motif_name,'\\('), function(x){
    strsplit(x[1],'/')[[1]][1]
  })%>%unlist()
  tf<-unique(tf)
  data.frame(tf=tf, dataset=x,num=1)
}))
saveRDS(mqtlsHomer,file.path(CONFIG$Result,'mqtls.homer.rds'))
mqtlsHomer<-readRDS(file.path(CONFIG$Result,'mqtls.homer.rds'))
table(mqtlsHomer$tf)%>%sort
table(mqtlsHomer$dataset)%>%sort

a<-dcast(mqtlsHomer,tf~dataset,value.var = "num",fill=0)
m<-as.matrix(a[,-1])
rownames(m)<-a$tf
m<-m[,match(sheets, colnames(m))]

saveImage("smr.homer.2.pdf",width = 6,height = 30)
Heatmap(m,
        col=colorRamp2(c(0, 1),c("#fffeee", "#c82423")),
        cluster_columns = FALSE
        )
dev.off()
saveImage("smr.homer.pdf",width = 5,height = 8)
Heatmap(m,
        col=colorRamp2(c(0, 1),c("#fffeee", "#c82423")),
        cluster_columns = FALSE
)
dev.off()

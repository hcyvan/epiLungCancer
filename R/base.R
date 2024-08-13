library(dplyr)
library(GenomicRanges)
library(readxl)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


# Configure --------------------------------------------------------------------
configDefaultPath <- './config.default.yaml'
configPath <- './config.yaml'
if (!file.exists(configPath)) {
  if(file.exists(configDefaultPath)) {
    configPath <- configDefaultPath
  } else {
    stop('please add configure file [config.yaml] to this project root directory')
  }
}

CONFIG <- yaml::read_yaml(configPath)
CONFIG['DataRaw']<-file.path(CONFIG$DataDir, 'raw')
CONFIG['DataInter']<-file.path(CONFIG$DataDir, 'intermediate')
CONFIG['DataResult']<-file.path(CONFIG$DataDir, 'result')
# Base Function ----------------------------------------------------------------
## IO --------------------------------------------------------------------------
getExt <- function(file_path) {
  file_name <- basename(file_path)
  parts <- strsplit(file_name, "\\.")[[1]]
  if (length(parts) > 1) {
    tolower(tail(parts, 1))
  } else {
    NULL
  }
}

loadData <- function(name, ext=NULL, header=FALSE, force.refresh=FALSE, no.cache=FALSE){
  if (is.null(ext)) {
    ext<-getExt(name)
  }
  if (is.null(ext)) {
    stop(sprintf("Unknown file type: %s", name))
  }
  ext.path <- name
  rds.path <- paste(ext.path,'rds',sep = '.')
  if(file.exists(rds.path) && !force.refresh){
    readRDS(rds.path)
  }else if(file.exists(ext.path)) {
    if(ext=='csv'){
      data<-read.csv(ext.path)
      if (!no.cache){
        saveRDS(data, rds.path)
      }
      data
    }else if(ext=='tsv'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE)
      if (!no.cache){
        saveRDS(data, rds.path)
      }
      data
    }else if(ext=='bed'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE, header=header)
      if (!no.cache){
        saveRDS(data, rds.path)
      }
      data
    } else {
      stop(sprintf("The file format <%s> is not supported", ext)); 
    }
  } else{
    stop(sprintf("File not exist: %s and %s", rds.path, ext.path))
  }
}

saveImage <- function(file,...){
  file.path=file.path(CONFIG$DataResult, file)
  if (endsWith(file, '.pdf')){
    pdf(file=file.path, ...)
  }
}

saveBed<-function(bed,outfile,...){
  if (!startsWith(colnames(bed)[1], '#')){
    colnames(bed)[1]<-paste0('#',colnames(bed)[1])
  }
  write.table(bed,outfile,sep = '\t',quote = FALSE,row.names = FALSE,...)
}
## Data ------------------------------------------------------------------------
removeNegativeOne <- function(m){
  m[rowSums(m[,4:ncol(m)]==-1)==0,]
}

#' print the confidence interval of an array
#'
#' @param data the target data array 
#' @param conf.level the confidence level
#'
#' @return
#' @export
#'
#' @examples
#' data <- c(2.3, 3.1, 2.8, 3.6, 3.2, 3.8, 3.0, 2.9, 3.4)
#' printConfidenceInterval(data, conf.level = 0.9) 
#' 
#' [1] "3.1222 (90% CI: 2.8437 - 3.4008)"
#' 
printConfidenceInterval<-function(data, conf.level=0.95){
  test<-t.test(data, conf.level = conf.level)
  v<-test$estimate
  conf.int<-test$conf.int
  conf<-attributes(conf.int)$conf.level*100
  sprintf("%.4f (%d%% CI: %.4f - %.4f)", v, conf,conf.int[1],conf.int[2])
}
## Genome Region ---------------------------------------------------------------
annoPeak<-function(peakGr) {
  peakAnno <- annotatePeak(peakGr, tssRegion = c(-5000, 5.000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  peakAnno@anno$symbol<-mapIds(org.Hs.eg.db, keys = peakAnno@anno$geneId, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  peakAnno@anno$ensembl<-mapIds(org.Hs.eg.db, keys = peakAnno@anno$geneId, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")
  anno<-peakAnno@anno
  anno[abs(anno$distanceToTSS)<100000,]
}
#' Calculate the intersection and difference of two GRanges
#'
#' @param gr0
#' @param gr1
#'
#' @return
#' @export
#'
#' @examples
#' 
grOverlap<-function(gr0,gr1){
  overlaps<-findOverlaps(gr0,gr1)
  i0<-unique(queryHits(overlaps))
  i1<-unique(subjectHits(overlaps))
  gr0i<-gr0[i0]
  gr0u<-gr0[-i0]
  gr1i<-gr1[i1]
  gr1u<-gr1[-i1]
  gr0stat<-unlist(
    list(
      gr0=length(gr0),
      gr0u=length(gr0u),
      gr0i=length(gr0i)
    )
  )
  gr1stat<-unlist(
    list(
      gr1=length(gr1),
      gr1u=length(gr1u),
      gr1i=length(gr1i)
    )
  )
  list(
    gr0stat=gr0stat,
    gr1stat=gr1stat,
    gr0=gr0,
    gr0u=gr0u,
    gr0i=gr0i,
    gr1=gr1,
    gr1u=gr1u,
    gr1i=gr1i
  )
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
# Project ----------------------------------------------------------------------
## Variable --------------------------------------------------------------------
FACTOR_LEVEL_GROUP<-c('CTL', 'LUAD', 'LUSC', 'LCC','SCLC')
COLOR_MAP_GROUP<-c('#2878b5','#c82423','#ffb15f','#fa66b3', '#925ee0')
names(COLOR_MAP_GROUP)<-FACTOR_LEVEL_GROUP
SHAPE_MAP_GROUP=c(16,15,17,10, 7)
names(SHAPE_MAP_GROUP)<-FACTOR_LEVEL_GROUP
HUMAN_CHROMSOME <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
  "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
  "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
  "chr22", "chrX", "chrY", "chrM"
)

## Function --------------------------------------------------------------------
g2color<-function(group,alpha=NULL){
  color<-COLOR_MAP_GROUP[match(group, names(COLOR_MAP_GROUP))]
  if(!is.null(alpha)){
    color<-alpha(color, alpha)
  }
  color
}
## Class -----------------------------------------------------------------------
### Sample ---------------------------------------------------------------------
Sample <- setRefClass(
  "sample",
  fields = list(table = "data.frame", list='list',color.map='data.frame'),
  methods = list(
    initialize = function(data) {
      data<-read_excel(file.path(CONFIG$DataRaw, 'SupplementaryData.xlsx'),sheet = 'sample')
      data<-data.frame(data)
      data$Group<-factor(data$Group, levels = FACTOR_LEVEL_GROUP)
      data<-arrange(data, Group)
      data$Color<-COLOR_MAP_GROUP[match(data$Group, names(COLOR_MAP_GROUP))]
      table<<-dplyr::select(data, SampleName, Group, Color)
      list<<-split(table,table$Group)
    },
    sample2group=function(samples){
      table$Group[match(samples,table$SampleName)]
    },
    selectFromBed=function(bed,samples=NULL){
      if (is.null(samples)){
        samples<-table$SampleName
      }
      bed[,c(1:3,match(samples, colnames(bed)))]
    },
    selectFromMatrixByGroup=function(data,group){
      df<-filter(table, Group==group)
      coli<-match(df$SampleName, colnames(data))
      data[,coli]
    },
    show = function() {
      print(table(table$Group))
    }
  )
)
## Class instance -------------------------------------------------------------
SAMPLE<-Sample()

# Plot function ----------------------------------------------------------------
plotBedListDensity<-function(bedList,class,target='tss',
                             lim=5000,
                             title="",
                             xlab="Distance",
                             color.map=COLOR_MAP_GROUP){
  genomicRegion<-readRDS(file.path(CONFIG$DataRaw,'genomicRegion.rds'))
  if (target=='tss'){
    target.gr<-genomicRegion$tss
  }else if (target=='cgi'){
    target.gr<-genomicRegion$cgIslands
  }

  dist<-lapply(bedList, function(x){
    bed<-x[[class]]
    distTo<-annoDistQueryToSubject(bed2gr(bed),target.gr)
    distTo
  })

  dens<-lapply(dist,function(x){
    d0<-x[abs(x$distanceToSubject)<lim,]
    density(d0$distanceToSubject)
  })
  xlim<-range(sapply(dens, "[", "x"))
  ylim<-range(sapply(dens, "[", "y"))
  plot(NA, main = title,
       xlab = xlab,
       ylab = "DMRs Density",
       xlim=xlim,
       ylim=ylim)
  for(i in 1:length(color.map)) {
    lines(dens[[names(color.map)[i]]], lwd = 2,col=color.map[i])
  }
  legend("topleft", legend=names(color.map), fill=color.map, bty = "n")
}

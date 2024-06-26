library(dplyr)
library(GenomicRanges)
library(readxl)


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
################################################################################
getExt <- function(file_path) {
  file_name <- basename(file_path)
  parts <- strsplit(file_name, "\\.")[[1]]
  if (length(parts) > 1) {
    tolower(tail(parts, 1))
  } else {
    NULL
  }
}

loadData <- function(name, ext=NULL, header=FALSE, force.refresh=FALSE){
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
      saveRDS(data, rds.path)
      data
    }else if(ext=='tsv'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE)
      saveRDS(data, rds.path)
      data
    }else if(ext=='bed'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE, header=header)
      saveRDS(data, rds.path)
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

removeNegativeOne <- function(m){
  m[rowSums(m[,4:ncol(m)]==-1)==0,]
}


feature2Bed<-function(feature) {
  bed<-do.call(rbind,lapply(feature, function(x){
    tmp<-strsplit(x,':')[[1]]
    se<-strsplit(tmp[2],'-')[[1]]
    c(tmp[1], se[1], se[2])
  }))
  bed<-data.frame(bed)
  bed[,2]<-as.numeric(bed[,2])
  bed[,3]<-as.numeric(bed[,3])
  colnames(bed)<-c('chrom','start','end')
  bed
}

feature2GRanges<-function(feature){
  bed2GRanges(feature2Bed(feature))
}

bed2Feature<-function(bed){
  paste(paste(bed[,1], bed[,2], sep=':'),bed[,3],sep='-')
}

bed2GRanges<-function(bed){
  gr<-GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2]+1, end =  bed[,3]))
  if (ncol(bed)>=4){
    .<-lapply(colnames(bed)[4:ncol(bed)],function(x){
      mcols(gr)[[x]]<<-bed[[x]]
    })
  }
  gr
}

GRanges2bed<-function(gr){
  data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
}
GRanges2Feature<-function(gr){
  bed<-data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
  bed2Feature(bed)
}
################################################################################
FACTOR_LEVEL_GROUP<-c('CTL', 'LUAD', 'LUSC', 'LCC','SCLC')
COLOR_MAP_GROUP<-c('#2878b5','#c82423','#ffb15f','#fa66b3', '#925ee0')
names(COLOR_MAP_GROUP)<-FACTOR_LEVEL_GROUP
SHAPE_MAP_GROUP=c(16,15,17,10, 7)
names(SHAPE_MAP_GROUP)<-FACTOR_LEVEL_GROUP



g2color<-function(group,alpha=NULL){
  color<-COLOR_MAP_GROUP[match(group, names(COLOR_MAP_GROUP))]
  if(!is.null(alpha)){
    color<-alpha(color, alpha)
  }
  color
}
#################################################################################
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
    show = function() {
      print(table(table$Group))
    }
  )
)

SAMPLE<-Sample()



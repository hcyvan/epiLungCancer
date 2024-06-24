#' Conversion between different genome coordinate formats
#'
#' Description: Conversion between different genome coordinate formats
#'   
#'
#' Details: There are 3 different genome coordinate formats is this document.
#'              1. gr: GRanges in GenomicRanges package. 1-base include.
#'              2. bed: data.frame, and thefirst three columns of the data.frame
#'                      are chrom, start, and end. 0-base exclude.
#'              3. feature: string like chr1:10-20. 0-base exclude.#'
#' @docType package
#' @name easyapi
NULL


#' Convert feature format to bed format
#' 
#' @param feature feature format data
#' @return bed format data
#' @examples
#' feature2bed
#' @export
feature2bed<-function(feature) {
  bed<-do.call(rbind,lapply(feature, function(x){
    tmp<-strsplit(x,':')[[1]]
    se<-strsplit(tmp[2],'-')[[1]]
    c(tmp[1], se[1], se[2])
  }))
  bed<-data.frame(bed)
  bed[,2]<-as.numeric(bed[,2])
  bed[,3]<-as.numeric(bed[,3])
  colnames(bed)<-c('chrom','start','end')
  return(bed)
}

#' Convert feature format to gr format
#'
#' @param feature feature format data
#'
#' @return gr format data
#' @export
#'
#' @examples
feature2gr<-function(feature){
  bed2gr(feature2bed(feature))
}

#' Convert bed format to feture format
#'
#' @param bed bed format data
#'
#' @return feature format data
#' @export
#'
#' @examples
bed2feature<-function(bed){
  paste(paste(bed[,1], bed[,2], sep=':'),bed[,3],sep='-')
}

#' Convert bed format to gr format
#'
#' @param bed bed format data
#'
#' @return gr format data
#' @export
#'
#' @examples
bed2gr<-function(bed){
  gr<-GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2]+1, end =  bed[,3]))
  if (ncol(bed)>=4){
    .<-lapply(colnames(bed)[4:ncol(bed)],function(x){
      mcols(gr)[[x]]<<-bed[[x]]
    })
  }
  gr
}

#' Convert gr format to bed format
#'
#' @param gr gr format data
#'
#' @return bed format data
#' @export
#'
#' @examples
gr2bed<-function(gr){
  data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
}

#' Convert gr format to feture format
#'
#' @param gr gr format data
#'
#' @return feature format data
#' @export
#'
#' @examples
gr2feature<-function(gr){
  bed<-data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
  bed2feature(bed)
}
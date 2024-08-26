library(dplyr)

#' Extract the results from the findMotifsGenome.pl script
#' 
#' Currently only knownResults.txt results are extracted. Filtering of results is supported
#'
#' @param result.dir the result directory of  findMotifsGenome.pl
#' @param p.value the maximum p-value of the results
#' @param q.value the maximum q-value of the results
#' @param foldchange the minimum fold change of the results 
#'
#' @return
#' @export
#'
#' @examples
#' getfindMotifsGenomeResults('/dmc/homer80/LUAD.hypo',p.value=0.05, foldchange=1.2)
#' 
getfindMotifsGenomeResults<-function(result.dir,p.value=NULL,q.value=NULL,foldchange=NULL){
  findMotifsGenome<-read.csv(file.path(result.dir, 'knownResults.txt'),sep='\t', check.names = FALSE)
  colnames(findMotifsGenome)<-c('motif_name','consensus','p','logp','q','target','target_percent','background','background_percent')
  findMotifsGenome$target_percent<-percent2numeric(findMotifsGenome$target_percent)
  findMotifsGenome$background_percent<-percent2numeric(findMotifsGenome$background_percent)
  findMotifsGenome$motif<-sapply(strsplit(findMotifsGenome$motif_name, '/'),function(x){
    x[1]
  })
  findMotifsGenome$tf<-sapply(strsplit(findMotifsGenome$motif_name, '/'),function(x){
    strsplit(x,"\\(")[[1]][1]
  })
  findMotifsGenome$fc <-findMotifsGenome$target_percent/findMotifsGenome$background_percent
  if (!is.null(p.value)){
    findMotifsGenome<-filter(findMotifsGenome, p<=p.value)
  }
  if (!is.null(q.value)){
    findMotifsGenome<-filter(findMotifsGenome, q<=q.value)
  }
  if (!is.null(foldchange)){
    findMotifsGenome<-filter(findMotifsGenome, fc>=foldchange)
  }
  return(findMotifsGenome)
}


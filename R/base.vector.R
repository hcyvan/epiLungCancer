source('R/base.R')
library(easyepi)
library(gtools)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggExtra)
library(ggplot2)
library(gridExtra)
library(ggforce)
library(patchwork)
library(karyoploteR)
library(rGREAT)

# Function ---------------------------------------------------------------------
interDmrVector<-function(dmr,vfile){
  suppressWarnings(
    {
      vector<-loadData(vfile,ext = 'bed', no.cache = TRUE)
      vector<-vector[,c(1,3,4,2,5:ncol(vector))]
      gr2<-bed2gr(vector)
      gr1<-bed2gr(dmr)
      overlaps<-findOverlaps(gr1,gr2)
      a1<-gr1[queryHits(overlaps)]
      a2<-gr2[subjectHits(overlaps)]
      a1<-unique(a1)
      a2<-unique(a2)
      out<-c(length(gr1), length(gr2), length(a1), length(a2))
      names(out)<-c('DMR', 'Vector', 'InterDMR','InterVector')
      out
    }
  )
}

readMV<-function(mvm.file){
  lines <- readLines(mvm.file)
  lines <- grep("^##", lines, invert = TRUE, value = TRUE)
  data <- read.csv(text = paste(lines, collapse = "\n"),check.names = FALSE,sep = '\t')
  return(data)
}
## class Motif ----------------------------------------------------------------
Motif <- setRefClass(
  "Motif",
  fields = list(count = "numeric", motifs='character'),
  methods = list(
    initialize = function(count) {
      count<<-count
      motifs<<-get_motif_array()
    },
    get_motif_array=function(type='CT'){
      if (type=='CT'){
        M<-'C'
        U<-'T'
      }else{
        M<-'1'
        U<-'0'
      }
      
      ms <- c()
      ms <- c(ms, paste(rep(M, count), collapse = ""))
      for (i in 1:count) {
        combos <- combinations(count, i)
        for (e in 1:nrow(combos)) {
          motif <- rep(M, count)
          indices <- combos[e, ] 
          motif[indices] <- U
          ms <- c(ms, paste(motif, collapse = ""))
        }
      }
      return(ms)
    },
    motif2vector = function(format='list') {
      m2v <- list()
      for (m in motifs) {
        vector <- sapply(strsplit(m, NULL)[[1]], function(x) ifelse(x == "C", 1, 0))
        m2v[[m]] <- vector
      }
      if(format=='list'){
        m2v
      }else{
        a<-do.call(rbind,m2v)%>%t
        rownames(a)<-NULL
        a
      }
    },
    show = function() {
      print(motifs)
    }
  )
)
## class Beta ------------------------------------------------------------------
Beta <- setRefClass(
  "Beta",
  fields = list(table = "data.frame",key='character'),
  methods = list(
    initialize = function(bfile) {
      beta<-read.csv(bfile,sep='\t')
      table<<-beta
      key<<-paste(paste(beta$X.chrom, beta$start,sep = ':'), beta$end, sep = '-')
    },
    getBeta=function(region){
      data<-table[match(region,key),]
      if (nrow(data)==1){
        unlist(data[,-1:-3])
      }else{
        data
      }
    }
  )
)
## class SMVC ------------------------------------------------------------------
SMVC <- setRefClass(
  "SMVC",
  fields = list(table = "data.frame",
                hypo = "data.frame",
                hyper = "data.frame",
                table.ori = "data.frame",
                reduceBed = 'list',
                reduceGr = 'list',
                filename2='character',
                filename='character'),
  methods = list(
    initialize = function(mvm.file,mvm.file2=NULL) {
      filename<<-mvm.file
      if (!is.null(mvm.file2)){
        filename2<<-mvm.file2
      }
      table.ori<<-readSMVC(mvm.file)
      if (!is.null(table.ori)){
        table<<-dplyr::select(table.ori, 'chrom', 'start', 'end', 'cpg', 'mvs', 'c_num','smvp','ssp','i0','i1','i0v0','i1v0','labels')
        table$class<<-ifelse(table$i0v0>table$i1v0, 'hypo','hyper')
        hypo<<-filter(table, class=='hypo')
        hyper<<-filter(table, class=='hyper')
        reduceGr<<-list(
          table=reduce(bed2gr(table),with.revma=TRUE),
          hypo=reduce(bed2gr(hypo),with.revma=TRUE),
          hyper=reduce(bed2gr(hyper),with.revma=TRUE)
        )
        reduceBed<<-lapply(reduceGr, function(x){
          out<-data.frame(gr2bed(x))
          out$revmap<-lapply(out$revmap,function(x){
            paste(x,collapse  = ',')
          })%>%unlist
          as.data.frame(out)
        })
      }
    },
    readSMVC = function(vfile){
      vector<-loadData(vfile,ext = 'bed')
      if (ncol(vector)==0){
        return(NULL) 
      }
      vector<-filter(vector, !V1%in%c('chrM','chrX','chrY'))
      names(vector)<-c('chrom','cpg','start','end','mvs','c_num','c_center','c_group_mvs_num','c_group_samples_num','c_group_samples','smvp','ssp','i0','i1','i0v0','i1v0','labels')
      vector[,c(1,3,4,2,5:ncol(vector))]
    },
    overlap = function(data, bed){
      #' gr0: regions in MVC
      #' gr0u: regions only in MVC
      #' gr0i: regions in MVC and DMR
      #' gr1: regions in DMR
      #' gr1u: regions only in DMR
      #' gr1i: regions in DMR and MVC
      gr0<-bed2gr(data)
      gr1<-bed2gr(bed)
      ov<-grOverlap(gr0,gr1)
      ov
    },
    .freeze=function(ouDir,dmr=NULL,class="hypo", do.reduce=FALSE){
      smvc.file<-sub("\\.smvc$", ".final.smvc", filename2)
      saveBed(table.ori[,c(1,4,2,3,5:ncol(table.ori))], smvc.file, col.names=FALSE)
      if (do.reduce){
        out.suffix<-'.reduce.bed'
        table.data<-reduceBed[[class]]
        table.data.table<-reduceBed$table
      }else{
        out.suffix<-'.bed'
        table.data<-.self[[class]]
        table.data.table<-.self$table
        # table.data<-table.data[,1:3]
        # table.data$feature<-bed2feature(table.data)
      }
      bed.file<-file.path(ouDir,sub("\\.smvc$", out.suffix, basename(filename2)))
      saveBed(table.data.table, bed.file,col.names=FALSE)
      if (!is.null(dmr)){
        ov<-overlap(table.data,dmr)
        for (x in c("gr0","gr0u","gr0i","gr1","gr1u","gr1i")){
          data<-ov[[x]]
          data<-gr2bed(data)
          bed.file<-file.path(ouDir,sub("\\.smvc$", paste0('.',class,'.',x,out.suffix), basename(filename2)))
          data$feature<-bed2feature(data)
          data<-dplyr::select(data, chrom, start,end,feature)
          saveBed(data, bed.file,col.names=FALSE)
          
        }
      }
    },
    freeze=function(outDir,dmr=NULL){
      .freeze(outDir, dmr$hypo, 'hypo', FALSE)
      .freeze(outDir, dmr$hypo, 'hypo', TRUE)
      .freeze(outDir, dmr$hyper, 'hyper', FALSE)
      .freeze(outDir, dmr$hyper, 'hyper', TRUE)
    },
    cpg2genome=function(cpgIdx){
      data<-filter(table, cpg==cpgIdx)
      paste0(data$chrom, ':',data$start, '-',data$end)
    },
    show = function() {
      print(table)
    }
  )
)
#smvc<-SMVC(file.path(CONFIG$DataInter,'vector/w4/LUAD_1.0_0.2.smvc'),file.path(CONFIG$DataInter,'vector/w4/LUAD.smvc'))
#smvc$freeze(file.path(CONFIG$DataInter,'vector/w4/bed'), dmrList$LUAD)

## class MVM ------------------------------------------------------------------
MVM <- setRefClass(
  "MVM",
  fields = list(table = "data.frame"),
  methods = list(
    initialize = function(mvm.file) {
      table<<-readMV(mvm.file) 
    },
    getBySample=function(sample, window=NULL){
      col<-table[sample]
      col<-unlist(col)
      names(col)<-table$cpg
      m<-mvs2matrix(col)
      if (!is.null(window)){
        motif<-Motif(4)
        colnames(m)<-motif$get_motif_array('01')
      }
      return(m)
    },
    getByCpG=function(cpg.idx){
      row<-table[table$cpg==cpg.idx,-1:-2]
      row<-unlist(row)
      return(mvs2matrix(row))
    },
    mvs2matrix=function(mvs){
      sapply(strsplit(mvs,'\\|'), function(x){
        as.numeric(x)
      })%>%t
    }
  )
)
## class MVH ------------------------------------------------------------------
MVH <- setRefClass(
  "MVH",
  fields = list(table = "data.frame"),
  methods = list(
    initialize = function(mv.file) {
      table<<-readMV(mv.file) 
    },
    getBySample=function(sample){
      col<-table[sample]
      return(unlist(col))
    },
    getByCpG=function(cpg.idx){
      row<-table[table$cpg==cpg.idx,-1:-2]
      return(unlist(row))
    },
    matrix=function(col.name.pattern=NULL){

      m<-table[,-1:-2]
      rownames(m)<-table$cpg
      m<-as.matrix(m)
      if (!is.null(col.name.pattern)){
        m<-m[,sapply(colnames(m), function(x){grepl(col.name.pattern,x)})]
      }
      return(m)
    }
  )
)
## plot function ---------------------------------------------------------------
plotMotifAnnoGgplot<-function(m=5, ori='v',reverse=FALSE){
  motif<-Motif(m)
  vcs<-motif$motif2vector()
  seqi<-1:length(vcs)
  sites<-lapply(1:length(vcs), function(y){
    xi<-sapply(1:length(vcs[[y]]),function(x){
      c(x, y,vcs[[y]][x])
    })%>%t
    xi
  })
  data<-data.frame(do.call(rbind,sites))
  data$motif<-motif$motifs[data$V2]
  data$motif<-factor(data$motif, levels = motif$motifs)
  data$r<-0.45
  data$C<-as.character(data$C)
  color_map<-c('black','white')
  names(color_map)<-c('1','0')
  if (ori=='v'){
    data<-data.frame(x0=data$V1,y0=data$V2, r=data$r, x=data$V1, y=data$motif, fill=data$C)
  }else if(ori=='h') {
    data<-data.frame(x0=data$V2,y0=data$V1, r=data$r, x=data$motif, y=data$V1, fill=data$C)
  }
  if (reverse){
    data$y0<-max(data$y0)-data$y0+1
  }
  ggplot(data, aes(x0 = x0, y0 = y0, r = r,x=x,y=y,fill=fill)) +
    geom_circle(color = "black")+
    scale_fill_manual(values = color_map) +
    theme_void()+
    # coord_fixed()+
    theme(legend.position = "none")
}

plotMotifAnnoHeatmap<-function(count, ori='h'){
  motif<-Motif(count)
  m = motif$motif2vector(format='mat')
  if (ori=='v'){
    m<-t(m)
    m<-m[nrow(m):1,]
    if (count==4){
      ws<-0.5
    }else{
      ws<-3
    }
    Heatmap(m, rect_gp = gpar(type = "none"),
            width = unit(2, "cm"),
            show_heatmap_legend = FALSE,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            show_column_names = FALSE,
            cell_fun = function(j, i, x, y, w, h, fill) {
              p = m[i, j]
              if (p==1){
                fill.color<-'black'
              }else{
                fill.color<-"white"
              }
              grid.circle(x,y, w*ws, gp = gpar(fill = fill.color, col = "black"))
            })
  }else{
    if (count==4){
      ws<-1.9
    }else{
      ws<-3
    }
    Heatmap(m, rect_gp = gpar(type = "none"),
            height = unit(2, "cm"),
            show_heatmap_legend = FALSE,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            show_column_names = FALSE,
            cell_fun = function(j, i, x, y, w, h, fill) {
              p = m[i, j]
              if (p==1){
                fill.color<-'black'
              }else{
                fill.color<-"white"
              }
              grid.circle(x,y, w*ws, gp = gpar(fill = fill.color, col = "black"))
            })
  }
}

plotMotifFrac<-function(data){
  ma<-sapply(split(SAMPLE$table,SAMPLE$table$Group), function(x){
    colSums(data[match(x$SampleName,rownames(data)),])
  })%>%t
  ma<-t(apply(ma, 2, function(x) x/sum(x)))
  ma[is.na(ma)]<-0
  ma<-ma*100
  
  p<-anno_barplot(ma, gp = gpar(fill = COLOR_MAP_GROUP), bar_width = 1, height = unit(2, "cm"))
  col_annotation = HeatmapAnnotation(p = p, show_annotation_name = FALSE)
  return(col_annotation)
}

plotWindowBase<-function(data,beta=NULL){
  row_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=SAMPLE$sample2group(rownames(data))),
    col = list(Stage =COLOR_MAP_GROUP),
    show_legend = FALSE,
    show_annotation_name =FALSE,
    which = 'row'
  )
  beta_annotation=NULL
  if (!is.null(beta)){
    p<-anno_points(beta,which='row',pch = 20,size = unit(0.8, "mm"),ylim = c(0,1))
    beta_annotation = HeatmapAnnotation(
      p = p,
      which = 'row',
      show_annotation_name = FALSE)
  }
  m<-as.matrix(data)
  # m<-log10(m+0.00001)
  h<-Heatmap(m,
             name='Log10(count)',
             col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
             left_annotation  = row_annotation,
             right_annotation = beta_annotation,
             # top_annotation =col_annotation,
             # height = unit(1, "npc"),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             show_column_names = FALSE
  )
  return(h)
}

plotWindowBase2<-function(data){
  m<-as.matrix(data)
  Heatmap(m,
          name='Count',
          col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
          # row_names_gp = gpar(fontsize = 5),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = TRUE
  )
}

plotWindow<-function(data,window=4,beta=NULL){
  annoMotif<-plotMotifAnnoHeatmap(window)
  annoFrac<-plotMotifFrac(data)
  h<-plotWindowBase(data,beta=beta)
  annoFrac%v%h%v%annoMotif
}
plotWindow2<-function(data,window=4){
  h<-plotWindowBase2(data)
  annoMotif<-plotMotifAnnoHeatmap(window)
  h%v%annoMotif
}


mvcReshape<-function(mvc){
  x<-0
  ret<-do.call(rbind,apply(mvc, 1, function(num){
    x<<-x+1
    y0<-1:16
    y<-ifelse(num!=0, y0,0)
    data.frame(x,num,y,y0)
  }))
  rownames(ret)<-NULL
  ret<-filter(ret, y!=0)
  ret
}

plotMvc<-function(mvc, window=4){
  data<-mvcReshape(mvc)
  motif<-Motif(window)
  data$motif<-motif$motifs[data$y]
  data$y<-factor(as.character(data$y),levels = as.character(unique(sort(data$y))))
  data$motif<-factor(data$motif, levels = rev(motif$motifs))
  data$num<-log(data$num)
  ggplot(data=data, aes(x=x, y=motif,color=num)) +
    geom_point()+
    theme_bw()+
    scale_color_gradient(low="#2878b5", high="#c82423")+
    labs(x=NULL)+
    # labs(x = "Methylation vector windows")+
    # scale_color_gradient2(midpoint=mid, low="blue", mid="white",high="red", space ="Lab" )
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
}

plotSmvcRegion<-function(mvc){
  p0<-plotMotifAnnoGgplot(4,reverse = TRUE)
  p1<-plotMvc(mvc,4)
  pp<-p0+p1+plot_layout(widths = c(1, 29))
  pp
}

plotSmvcHeatmap<-function(mvc,window=4,group=NULL,mark=NULL,beta=NULL){
  mvc<-log(mvc+0.00001)
  m<-t(mvc)
  #m<-m[nrow(m):1,]
  row_annotation<-NULL
  if (!is.null(group)){
    row_annotation <-HeatmapAnnotation(
      df=data.frame(Stage=rep(group, nrow(m))),
      col = list(Stage =COLOR_MAP_GROUP),
      show_annotation_name =FALSE,
      show_legend = FALSE,
      which = 'row'
    )
  }
  fonts<-rep(0, ncol(m))
  if (!is.null(mark)){
    fonts[match(mark, colnames(m))]<-5
  }
  beta_annotation=NULL
  if (!is.null(beta)){
    p<-anno_points(beta,which='col',pch = 20,size = unit(0.2, "mm"),ylim=c(0,1))
    beta_annotation = HeatmapAnnotation(
      p = p,
      which = 'col',
      show_annotation_name = FALSE)
  }
  Heatmap(m,
          show_heatmap_legend = TRUE,
          right_annotation  = row_annotation,
          top_annotation = beta_annotation,
          column_names_gp = gpar(fontsize = fonts),
          col=colorRamp2(c(0,max(m)), c("#ffffff", "#c82423")),
          show_column_dend =FALSE,
          cluster_rows = FALSE,
          cluster_columns = TRUE,
          clustering_method_columns = 'ward.D2',
          row_names_side = "left",
          show_row_names = TRUE,
          show_column_names = TRUE
  )
}

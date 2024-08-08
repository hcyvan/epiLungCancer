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
                filename='character'),
  methods = list(
    initialize = function(mvm.file) {
      filename<<-mvm.file
      table.ori<<-readSMVC(mvm.file)
      if (!is.null(table.ori)){
        table<<-dplyr::select(table.ori, 'chrom', 'start', 'end', 'cpg', 'mvs', 'c_num','smvp','ssp','i0','i1','i0v0','i1v0','labels')
        table$class<<-ifelse(table$i0v0>table$i1v0, 'hypo','hyper')
        hypo<<-filter(table, class=='hypo')
        hyper<<-filter(table, class=='hyper')
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
    toBed = function(){
      bed<-sub("\\.smvc$", ".bed", filename)
      saveBed(table, bed,col.names=FALSE)
      smvc<-sub("\\.smvc$", ".smvc", filename)
      saveBed(table.ori[,c(1,4,2,3,5:ncol(vector))], smvc,col.names=FALSE)
    },
    overlap = function(bed,class=NULL){
      #' gr0: regions in MVC
      #' gr0u: regions only in MVC
      #' gr0i: regions in MVC and DMR
      #' gr1: regions in DMR
      #' gr1u: regions only in DMR
      #' gr1i: regions in DMR and MVC
      data<-table
      if (!is.null(class)){
        if(class=='hyper'){
          data<-hyper
        }else if (class=='hypo'){
          data<-hypo
        }
      }
      gr0<-bed2gr(data)
      gr1<-bed2gr(bed)
      ov<-grOverlap(gr0,gr1)
      ov
    },
    saveOverlap=function(bed,class){
      ov<-overlap(bed, class)
      for (x in c("gr0","gr0u","gr0i","gr1","gr1u","gr1i")){
        data<-ov[[x]]
        data<-gr2bed(data)
        out.suffix<-paste0('.',class,'.',x,'.bed')
        out<-sub("\\.smvc$", out.suffix, filename)
        saveBed(data, out,col.names=FALSE)
      }
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

#smvc<-SMVC(file.path(CONFIG$DataInter,'vector/w4/LUAD_1.0_0.4.smvc'))
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
    matrix=function(mvs){
      m<-table[,-1:-2]
      rownames(m)<-table$cpg
      return(as.matrix(m))
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
  m<-log10(m+0.00001)
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

plotWindow<-function(data,window=5,beta=NULL){
  annoMotif<-plotMotifAnnoHeatmap(window)
  annoFrac<-plotMotifFrac(data)
  h<-plotWindowBase(data,beta=beta)
  annoFrac%v%h%v%annoMotif
}
plotWindow2<-function(data,window=5){
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
    scale_color_gradient(low="blue", high="red")+
    labs(x=NULL)+
    # labs(x = "Methylation vector windows")+
    # scale_color_gradient2(midpoint=mid, low="blue", mid="white",high="red", space ="Lab" )
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
}
plotMvc(mvm$getBySample('LUAD'))



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


# SMVC -------------------------------------------------------------------
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
groups<-names(COLOR_MAP_GROUP)[-1]
files<-sapply(groups, function(x){
  file.path(CONFIG$DataInter,'vector',window,paste0(x,'_1.0_0.4.smvc'))
})
smvcs<-lapply(files, function(x){
  SMVC(x)
})
## Data cleaning and statistics
### SMVC to bed ----------------------------------------------------------------
._<-lapply(smvcs, function(smvc){
  smvc$toBed()
})

### Intersect with DMR ----------------------------------------------------------------
do.call(rbind,lapply(groups, function(x){
  dmr<-dmrList[[x]]
  smvc<-smvcs[[x]]
  ov1<-smvc$overlap(dmr$hypo,class='hypo')
  ov2<-smvc$overlap(dmr$hyper,class='hyper')
  data<-rbind(c(ov1$gr0stat, ov1$gr1stat),c(ov2$gr0stat, ov2$gr1stat))
  colnames(data)<-c('mvc', 'mvcUniq', 'mvcInter','dmr','dmrUniq','dmrInter')
  data<-data.frame(group=x,data,class=c('hypo','hyper'))
  data
}))
### Save MVC and DMR -------------------------------------------------------------
._<-lapply(groups, function(x){
  dmr<-dmrList[[x]]
  smvc<-smvcs[[x]]
  smvc$saveOverlap(dmr$hypo, 'hypo')
  smvc$saveOverlap(dmr$hyper, 'hyper')
})
## Visualization -----------------------------------------------------------------
### LUAD ------------------------------------------------------------------------
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
bvalue<-beta$getBeta(smvc$cpg2genome(cpg))
data<-mvm$getByCpG(cpg)
saveImage("smvc.SCLC.hypo.4505825.pdf",width = 5,height = 6)
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
mark<-c(4505825, 16373819)
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/mvm/SCLC.group.mvm'))
p0<-plotSmvcHeatmap(mvm$getBySample('CTL',4),4,'CTL',mark,beta$table$CTL)
p1<-plotSmvcHeatmap(mvm$getBySample('LUAD',4),4,'LUAD',mark,beta$table$LUAD)
p2<-plotSmvcHeatmap(mvm$getBySample('LUSC',4),4,'LUSC',mark,beta$table$LUSQ)
p3<-plotSmvcHeatmap(mvm$getBySample('LCC',4),4,'LCC',mark,beta$table$LCC)
p4<-plotSmvcHeatmap(mvm$getBySample('SCLC',4),4,'SCLC',mark,beta$table$SCLC)
saveImage("smvc.SCLC.pdf",width =12,height = 6)
p0%v%p1%v%p2%v%p3%v%p4
dev.off()

#### Figure 5xxx smvc region ----------------------------------------------------
#mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.4505825.chr1_4505805_4505845.group.mvm'))
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.16373819.chr10_16373799_16373839.group.mvm')) # OK
#mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/region/SCLC.4505825.chr1_4505775_4505875.group.mvm')) # ok

pp0<-plotSmvcRegion(mvm$getBySample('CTL'))
pp1<-plotSmvcRegion(mvm$getBySample('LUAD'))
pp2<-plotSmvcRegion(mvm$getBySample('LUSC'))
pp3<-plotSmvcRegion(mvm$getBySample('LCC'))
pp4<-plotSmvcRegion(mvm$getBySample('SCLC'))
pp0/pp1/pp2/pp3/pp4

# HOMER analysis ---------------------------------------------------------------
fc<-1.3
q<-0.01

LUAD.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUAD_1.0_0.4.hypo.gr0'),q, fc)
LUAD.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUAD_1.0_0.4.hyper.gr0'),q, fc)
LUSC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUSC_1.0_0.4.hypo.gr0'),q, fc)
LUSC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LUSC_1.0_0.4.hyper.gr0'),q, fc)
LCC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LCC_1.0_0.4.hypo.gr0'),q, fc)
LCC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/LCC_1.0_0.4.hyper.gr0'),q, fc)
SCLC.hypo<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/SCLC_1.0_0.4.hypo.gr0'),q, fc)
SCLC.hyper<-getfindMotifsGenomeResults(file.path(CONFIG$DataInter,'vector/w4/homer/SCLC_1.0_0.4.hyper.gr0'),q, fc)

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
write.csv(output,file.path(CONFIG$DataInter,'vector/w4/homer','smvc.homer.csv'),row.names = FALSE)









# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx------
### CpG ------------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w4/mvm/LUAD.sample.mvm'))
data<-mvm$getByCpG(10958214) # methy level不明显
#saveImage("mv.window.pdf",width = 4.6,height = 6)
plotWindow(data,4)
#dev.off()


## GSE186458 -------------------------------------------------------------------
### CpG ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/LUAD.sample.mvm'))
data<-mvm$getByCpG(497192)
plotWindow2(data)

### Sample ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/LUAD.sample.mvm'))
data<-mvm$getBySample('Lung-Alveolar-Epithelial-Z000000T1')
plotWindow2(data)






















# All other human cell type ----------------------------------------------------


mvh.file<-file.path(CONFIG$DataInter,'vector/GSE79279/GSE79277/smvh/LUAD.smvh')
mvh.file<-file.path(CONFIG$DataInter,'vector/GSE79279/GSE79277/smvh/SCLC.smvh')
mvh.file<-file.path(CONFIG$DataInter,'vector/GSE79279/GSE79277/smvh/LCC.smvh')
mvh.file<-file.path(CONFIG$DataInter,'vector/GSE79279/GSE79277/smvh/LUSC.smvh')



# mvh.file<-file.path(CONFIG$DataInter,'vector/GSE186458/w5/mvh','LUAD.smvh')
# # mvh.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','LUSC.smvh')
# # mvh.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','LCC.smvh')
# mvh.file<-file.path(CONFIG$DataInter,'vector/GSE186458/w5/mvh','SCLC.smvh')
# mvh.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','LUAD.smvh')

# mvh.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','SCLC_1.4.smvh')
# mvh.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','LUAD_1.4.smvh')

mvh<-MVH(mvh.file)

m<-mvh$matrix()
m<-m[rowSums(m!=0)!=0,]

m<-log(m+0.0001)%>%t
# m<-mvh$matrix()%>%t%>%log
#saveImage("aaa.sclc.pdf",width = 49,height = 30)
#plotWindowBase2(m)
#dev.off()

Heatmap(m,
        name='Count',
        col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
        # row_names_gp = gpar(fontsize = 5),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE
)



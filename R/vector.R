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

# Function -----------------------------------------------------------------------------------------------------
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

interDmrVector2<-function(dmr, vector){
  suppressWarnings(
    {
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

Motif <- setRefClass(
  "Motif",
  fields = list(count = "numeric", motifs='character'),
  methods = list(
    initialize = function(count) {
      count<<-count
      motifs<<-get_motif_array()
    },
    get_motif_array=function(){
      ms <- c()
      ms <- c(ms, paste(rep("C", count), collapse = ""))
      for (i in 1:count) {
        combos <- combinations(count, i)
        for (e in 1:nrow(combos)) {
          motif <- rep("C", count)
          indices <- combos[e, ] 
          motif[indices] <- "T"
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

readMVM<-function(mvm.file){
  lines <- readLines(mvm.file)
  lines <- grep("^##", lines, invert = TRUE, value = TRUE)
  data <- read.csv(text = paste(lines, collapse = "\n"),check.names = FALSE,sep = '\t')
  return(data)
}

MVM <- setRefClass(
  "MVM",
  fields = list(table = "data.frame"),
  methods = list(
    initialize = function(mvm.file) {
      table<<-readMVM(mvm.file) 
    },
    getBySample=function(sample){
      col<-table[sample]
      col<-unlist(col)
      names(col)<-table$cpg
      return(mvs2matrix(col))
    },
    getByCpG=function(cpg.idx){
      row<-table[mvm$table$cpg==cpg.idx,-1:-2]
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


# Figure S3 ------------------------------

plotFrac<-function(group="LUAD", window='w4'){
  vfile<-file.path(CONFIG$DataInter,'vector',window,paste0(group,'_0.80_0.1.mvc'))
  vector<-loadData(vfile,ext = 'bed')
  vector<-vector[,c(1,3,4,2,5:ncol(vector))]
  vector<-(filter(vector, V11>=0.9, V12>0.3))
  data<-data.frame(
    x=vector$V11*100,
    y=vector$V12*100
  )
  p_scatter <- ggplot(data, aes(x = x, y = y)) +
    geom_point(size=0.3,alpha = 0.6) +
    # geom_density_2d_filled()+
    # xlab("MVs Frac")+
    # ylab("Sample Frac")+
    #ggtitle(paste0(group," ", window," ", nrow(data)))+
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

p2.3<-plotFrac('LUSC', 'w3')
p2.4<-plotFrac('LUSC', 'w4')
p2.5<-plotFrac('LUSC', 'w5')
p2.6<-plotFrac('LUSC', 'w6')
p2.7<-plotFrac('LUSC', 'w7')

p3.3<-plotFrac('LCC', 'w3')
p3.4<-plotFrac('LCC', 'w4')
p3.5<-plotFrac('LCC', 'w5')
p3.6<-plotFrac('LCC', 'w6')
p3.7<-plotFrac('LCC', 'w7')

p4.3<-plotFrac('SCLC', 'w3')
p4.4<-plotFrac('SCLC', 'w4')
p4.5<-plotFrac('SCLC', 'w5')
p4.6<-plotFrac('SCLC', 'w6')
p4.7<-plotFrac('SCLC', 'w7')

pp<-grid.arrange(
  p1.3,p1.4,p1.5,p1.6,p1.7,
  p2.3,p2.4,p2.5,p2.6,p2.7,
  p3.3,p3.4,p3.5,p3.6,p3.7,
  p4.3,p4.4,p4.5,p4.6,p4.7,
  ncol=5)
ggsave(file.path(CONFIG$DataResult, 'mv.dist.png'), plot = pp, width = 8, height = 8, dpi = 300)



# DMR vs SMV ------------------------------
dmrList<-readRDS(file.path(CONFIG$DataInter, 'dmc','p80','one2rest80.dmr.list.rds'))
loadMVC<-function(window, group){
  vfile<-file.path(CONFIG$DataInter,'vector',window,paste0(group,'_0.80_0.1.mvc'))
  vector<-loadData(vfile,ext = 'bed')
  vector<-vector[,c(1,3,4,2,5:ncol(vector))]
  vector<-filter(vector, !V1%in%c('chrM','chrX','chrY'))
}


window<-'w5'
sample.frac<-0.4
v0<-loadMVC(window,'CTL')
v1<-loadMVC(window,'LUAD')
v2<-loadMVC(window,'LUSC')
v3<-loadMVC(window,'LCC')
v4<-loadMVC(window,'SCLC')

mvsp<-list(
  CTL=filter(v0, V11>=1, V12>=sample.frac),
  LUAD=filter(v1, V11>=1, V12>=sample.frac),
  LUSC=filter(v2, V11>=1, V12>=sample.frac),
  LCC=filter(v3, V11>=1, V12>=sample.frac),
  SCLC=filter(v4, V11>=1, V12>=sample.frac)
)

interList<-list(
  hypoCTL=interDmrVector2(dmrList$CTL$hypo,mvsp$CTL),
  hyperCTL=interDmrVector2(dmrList$CTL$hyper,mvsp$CTL),
  hypoLUAD=interDmrVector2(dmrList$LUAD$hypo,mvsp$LUAD),
  hyperLUAD=interDmrVector2(dmrList$LUAD$hyper,mvsp$LUAD),
  hypoLUSC=interDmrVector2(dmrList$LUSC$hypo,mvsp$LUSC),
  hyperLUSC=interDmrVector2(dmrList$LUSC$hyper,mvsp$LUSC),
  hypoLCC=interDmrVector2(dmrList$LCC$hypo,mvsp$LCC),
  hyperLCC=interDmrVector2(dmrList$LCC$hyper,mvsp$LCC),
  hypoSCLC=interDmrVector2(dmrList$SCLC$hypo,mvsp$SCLC),
  hyperSCLC=interDmrVector2(dmrList$SCLC$hyper,mvsp$SCLC)
)
do.call(rbind,interList)



dmr1<-dmrList$LUAD$hypo
dmr2<-dmrList$LUAD$hyper



getOnco<-function(vector,dmr){
  hypo<-dmr$hypo
  hyper<-dmr$hyper
  gr0<-bed2gr(vector)
  gr1<-bed2gr(hypo)
  gr2<-bed2gr(hyper)
  c1<-gr0%>%length()
  c2<-setdiff(gr0,gr1)%>%length()
  c3<-setdiff(setdiff(gr0,gr1), gr2)%>%length()
  uni<-setdiff(setdiff(gr0,gr1), gr2)
  gr2bed(uni)
}


ocmv<-list(
  CTL=getOnco(mvsp$CTL, dmrList$CTL),
  LUAD=getOnco(mvsp$LUAD, dmrList$LUAD),
  LUSC=getOnco(mvsp$LUSC, dmrList$LUSC),
  LCC=getOnco(mvsp$LCC, dmrList$LCC),
  SCLC=getOnco(mvsp$SCLC, dmrList$SCLC)
)








for (group in names(mvsp)){
  vector<-mvsp[[group]]
  vector$V2<-as.character(vector$V2)
  vector$V3<-as.character(vector$V3)
  vector$V4<-as.character(vector$V4)
  mvc<-vector[,c(1,4,2,3,5:ncol(vector))]
  saveBed(mvc, file.path(CONFIG$DataInter,'vector',window,paste0(group,'.mvc')),col.names=FALSE)
  bed<-vector[,c(1,2,3,4,5)]
  saveBed(bed, file.path(CONFIG$DataInter,'vector',window,paste0(group,'.bed')),col.names=FALSE)
  saveBed(ocmv[[group]], file.path(CONFIG$DataInter,'vector',window,paste0(group,'.onco.bed')),col.names=FALSE)
}

# Heatmap of SMV  ---------------------------------------------------------------------------------------------
getHeatmapMotif<-function(count){
  motif<-Motif(count)
  m = motif$motif2vector(format='mat')
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
            grid.circle(x,y,
                        w*3,
                        gp = gpar(fill = fill.color, col = "black"))
          })
}

getHeatmapMotifFrac<-function(data){
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

getHeatmapInWindow<-function(data){
  row_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=SAMPLE$sample2group(rownames(data))),
    col = list(Stage =COLOR_MAP_GROUP),
    show_legend = FALSE,
    show_annotation_name =FALSE,
    which = 'row'
  )
  m<-as.matrix(data)
  h<-Heatmap(m,
             name='Count',
             col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
             left_annotation  = row_annotation,
             # top_annotation =col_annotation,
             # height = unit(1, "npc"),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             show_column_names = FALSE
  )
  return(h)
}

getHeatmapInWindow2<-function(data){
  m<-as.matrix(data)
  h<-Heatmap(m,
             name='Count',
             col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
             # row_names_gp = gpar(fontsize = 5),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = TRUE,
             show_column_names = FALSE
  )
}

plotMVwindow<-function(data){
  h<-getHeatmapInWindow(data)
  annoMotif<-getHeatmapMotif(5)
  annoFrac<-getHeatmapMotifFrac(data)
  annoFrac%v%h%v%annoMotif
}
plotMVwindow2<-function(data){
  h<-getHeatmapInWindow2(data)
  annoMotif<-getHeatmapMotif(5)
  h%v%annoMotif
}
## Lung cancer -----------------------------------------------------------------
### CpG ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/w5/mvm/LUAD.sample.mvm'))
data<-mvm$getByCpG(497192)
#saveImage("mv.window.pdf",width = 4.6,height = 6)
plotMVwindow(data)
#dev.off()

## GSE186458 -------------------------------------------------------------------
### CpG ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/LUAD.sample.mvm'))
data<-mvm$getByCpG(497192)
plotMVwindow2(data)

### Sample ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/LUAD.sample.mvm'))
data<-mvm$getBySample('Lung-Alveolar-Epithelial-Z000000T1')
plotMVwindow2(data)



# SMV in Genome ----------------------------------------------------------------
getMVsAnno<-function(m=5){
  motif<-Motif(m)
  vcs<-motif$motif2vector()
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
  p0<-ggplot(data, aes(x0 = V1, y0 = V2, r = r,x=V1,y=motif,fill=C)) +
    geom_circle(color = "black")+
    scale_fill_manual(values = color_map) +
    theme_void()+
    # coord_fixed()+
    theme(legend.position = "none")
  p0
}

getMVs<-function(mvc,n=5){
  motif<-Motif(5)
  sites<-lapply(1:nrow(mvc), function(i){
    mvs<-mvc$V5[i]
    num<-as.numeric(strsplit(mvs,'\\|')[[1]])
    sapply(1:length(num), function(j){
      y<-ifelse(num[j]==0, 0, j)
      c(i,y,num[j],j)
    })%>%t
  })
  data<-data.frame(do.call(rbind,sites))
  data<-filter(data, X2!=0)
  
  data$motif<-motif$motifs[data$X2]
  data$X2<-factor(as.character(data$X2),levels = as.character(unique(sort(data$X2))))
  data$motif<-factor(data$motif, levels = motif$motifs)
  
  data$X3<-log(data$X3)
  
  p1<-ggplot(data=data, aes(x=X1, y=motif,color=X3)) +
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
getBeta<-function(group){
  data<-loadData(bfile, header = TRUE,no.cache=TRUE)
  beta<-data[[group]]
  data<-data.frame(x=1:nrow(data),y=beta)
  p2<-ggplot(data=data, aes(x=x, y=y)) +
    # geom_bar(stat = "identity")+
    geom_point()+
    labs(x=NULL)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  p2
}


plotMVC<-function(vfile){
  mvc<-loadData(vfile,ext = 'bed', no.cache = TRUE)
  #idx<-sample(1:nrow(mvc),100)
  mvc<-mvc[1:100,]
  p0<-getMVsAnno(5)
  p1<-getMVs(mvc,5)
  pp<-p0+p1+plot_layout(widths = c(1, 29))
  pp
}

plotMVC2<-function(vfile, group){
  mvc<-loadData(vfile,ext = 'bed', no.cache = TRUE)
  p0<-getMVsAnno(5)
  p1<-getMVs(mvc,5)
  p2<-getBeta(group)
  pp<-(p2/p1)
  pp
}


#vfile<-file.path(CONFIG$DataInter,'vector','w5','LUAD.mvc')

pp0<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','mvm','LUAD/CTL.mvc'))
pp1<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','mvm','LUAD/LUAD.mvc'))
pp2<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','mvm','LUAD/LUSC.mvc'))
pp3<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','mvm','LUAD/LCC.mvc'))
pp4<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','mvm','LUAD/SCLC.mvc'))


pp0/pp1/pp2/pp3/pp4


pp0<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','group','SCLC/CTL.mvc'))
pp1<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','group','SCLC/LUAD.mvc'))
pp2<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','group','SCLC/LUSC.mvc'))
pp3<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','group','SCLC/LCC.mvc'))
pp4<-plotMVC(file.path(CONFIG$DataInter,'vector','w5','group','SCLC/SCLC.mvc'))
pp0/pp1/pp2/pp3/pp4

vfile<-file.path(CONFIG$DataInter,'vector','w5','group','LUAD/CTL.mvc')
mvc<-loadData(vfile,ext = 'bed', no.cache = TRUE)




# All other human cell type ----------------------------------------------------


mvm.file<-file.path(CONFIG$DataInter,'vector/w5/mvh','aa.mvh')

data<-readMVM(mvm.file)
rownames(data)<-data$cpg
m<-data[,3:ncol(data)]

print(dim(m))
m<-t(m)
m<-log(m)
# m<-m[100:200,]
# saveImage('aa.pdf',width = 16,height = 40)
Heatmap(m,
        name='Count',
        col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
        # left_annotation  = row_annotation,
        # top_annotation =col_annotation,
        # height = unit(1, "npc"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE
)
# dev.off()


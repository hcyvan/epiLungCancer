source('R/base.R')

library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
library(ComplexHeatmap)
library(dendextend)
library(ape)
library(FactoMineR)
library(factoextra)
library(Rtsne)

#----------------------------------------------------------------------------------------------------------------------
# Figure 1B: t-SNE analysis of all lung cancer samples
#----------------------------------------------------------------------------------------------------------------------
methy10x<-readRDS(file.path(CONFIG$DataRaw, 'merge.d10.one.bed.rds'))
set.seed(100)
idx<-sample(1:nrow(methy10x),200000)
m<-methy10x[idx, 4:ncol(methy10x)]
m<-na.omit(m)
m<-t(m)
tsne_result <- Rtsne(m, dims = 2, perplexity = 10, verbose = TRUE,max_iter = 1000)
tsne_df <- data.frame(X = tsne_result$Y[,1], Y = tsne_result$Y[,2], Species = SAMPLE$table$Group)
saveImage("methylation.level.tsne.pdf",width = 6,height = 5)
ggplot(tsne_df, aes(x = X, y = Y, color = Species, shape=Species)) +
  scale_color_manual(values = COLOR_MAP_GROUP)+
  geom_point(size = 2.5) +
  scale_shape_manual(values=SHAPE_MAP_GROUP)+
  theme_classic()+
  xlab('t-SNE 1')+
  ylab('t-SNE 2')
dev.off()
#==============================================================================================================
# p.s. PAC
#==============================================================================================================
methy10x<-readRDS(file.path(CONFIG$DataRaw, 'merge.d10.one.bed.rds'))
set.seed(100)
idx<-sample(1:nrow(methy10x),100000)
m<-methy10x[idx, 4:ncol(methy10x)]
res.pca <- PCA(t(((as.matrix(m)))), graph = FALSE)
saveImage("methylation.level.pca.pdf",width = 6,height = 5)
fviz_pca_ind(res.pca,repel = TRUE,label="none",col.ind=SAMPLE$table$Group,palette=COLOR_MAP_GROUP,pointsize =2.2)+
  theme_classic()+
  theme(plot.title = element_blank())
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# The confidence interval of bisufite conversion ratio and CpGs depth
#----------------------------------------------------------------------------------------------------------------------
data<-read_excel(file.path(CONFIG$DataRaw, 'SupplementaryData.xlsx'),sheet = 'qc')
printConfidenceInterval(bisulfite_conversion_ratio)
printConfidenceInterval(cg_depth)
#----------------------------------------------------------------------------------------------------------------------
# Figure 1C: The Average DNA methylation level of samples in each group
#----------------------------------------------------------------------------------------------------------------------
data<-read_excel(file.path(CONFIG$DataRaw, 'SupplementaryData.xlsx'),sheet = 'qc')

samplesMatch<-data[match(SAMPLE$table$SampleName,data$SampleName),]
data<-data.frame(
  group=SAMPLE$table$Group,
  ratio=samplesMatch$MCALL_MeanRatioCG_3X
)

saveImage("methylation.level.mean.pdf",width = 3.5,height = 2.5)
ggplot(data=data,aes(x=group,y=ratio,fill=group))+
  scale_fill_manual(values=COLOR_MAP_GROUP) +
  # geom_violin(trim=FALSE) +
  geom_violin( colour = 'NA',trim=TRUE)+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="#ffffff",size=0.2)+
  theme_classic()+
  stat_compare_means( comparisons = list(c('CTL','LUAD'),c('CTL','LUSC'),c('CTL','SCLC'),c('CTL','LCC')),
                      label = 'p.signif', method = "t.test")+
  labs(x='',y='Mean Methylation Level')+
  theme(legend.position="none",
        axis.title.x = element_text(size=0),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size = 14,colour="black"))+
  guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1D: The methylation level density distribution of CTL.vs.LUAD, CTL.vs.LUAC, CTL.vs.LCC and CTL.vs.SCLC
#----------------------------------------------------------------------------------------------------------------------
ratio <- loadData(file.path(CONFIG$DataRaw, 'merge.group.3x.all.bed'),header = TRUE)
ratio<-removeNegativeOne(ratio)
drawDensity<-function(ratio1,ratio2,s1,s2){
  d1 <- density(ratio1)
  d2 <- density(ratio2)
  dens <- list(a=d1,b=d2)
  plot(NA, xlim=range(sapply(dens, "[", "x")),
       ylim=range(sapply(dens, "[", "y")),xlab="",ylab="Density",cex.lab=1.5, cex.axis=1.5)
  polygon(d1, col=s1,border = NA)
  polygon(d2, col=s2,border = NA)
}
saveImage("methylation.level.density.pdf",width = 6,height = 4)
par(mfrow=c(4, 1))
par(mai=c(0, 0.5, 0, 0.5))
drawDensity(ratio$CTL,ratio$LUAD,g2color('CTL',0.7),g2color('LUAD',0.7))
drawDensity(ratio$CTL,ratio$LUSQ,g2color('CTL',0.7),g2color('LUSC',0.7))
drawDensity(ratio$CTL,ratio$LCC,g2color('CTL',0.7),g2color('LCC',0.7))
drawDensity(ratio$CTL,ratio$SCLC,g2color('CTL',0.7),g2color('SCLC',0.7))
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1E. Methylation level of CpGs within 5,000 bp upstream and downstream relative to CGI and TSS
#----------------------------------------------------------------------------------------------------------------------
#==============================================================================================================
# Part 1
#
# The signal.cgi.matrix.bed and signal.tss.matrix.bed are generate by `methytools signal`
#==============================================================================================================
smooth2 <- function(hw=51) {
  function(arr){
    for(i in 1:length(arr)){
      left<-i-hw
      if(left < 1){
        left<-1
      }
      right<-i+hw
      if(right > length(arr)){
        right<-length(arr)
      }
      
      arr[i]<-mean(arr[left:right])
    }
    arr
  }
}

plot.group.methy.profile<- function(UP,DWON,group,signal.file,ylab="Methylation Level",xlab="Distance to TSS (bp)",cex=1.5,horiz=FALSE, draw.legend=TRUE){
  x<-seq(UP, DWON)
  signalMatrixGroup<-loadData(signal.file,header = TRUE)
  mmData<-signalMatrixGroup[,-1]
  mm<-mmData[(10001+UP):(10001+DWON),]
  # mm<-t(scale(t(mm)))
  mm<-as.data.frame(apply(mm,2,smooth2(51)))
  mm<-as.data.frame(mm)
  par(mar = c(5,5,1,1))
  if (draw.legend){
    plot(NA, xlim=c(UP, DWON), xlab=xlab,ylab=ylab,ylim=c(min(mm),max(mm)), bty='n')
    box()
  }else{
    plot(NA, xlim=c(UP, DWON),cex.lab=1.5, cex.axis=1.5,ylim=c(min(mm),max(mm)), xaxt='n',yaxt='n',xlab = "",ylab = "")
  }
  lines(x,mm[["CTL"]],col=g2color('CTL',alpha = 0.7))
  lines(x,mm[[group]],col=g2color(group,alpha = 0.7))
}

saveImage("methylation.level.profile.cgi.pdf",width = 8,height = 2)
par(mfrow=c(1, 4))
plot.group.methy.profile(-4000, 4000,'LUAD',file.path(CONFIG$DataInter, "signal.cgi.matrix.bed"),xlab="Distance to CGI Center (bp)")
plot.group.methy.profile(-4000, 4000,'LUSC',file.path(CONFIG$DataInter, "signal.cgi.matrix.bed"),xlab="Distance to CGI Center (bp)")
plot.group.methy.profile(-4000, 4000,'LCC',file.path(CONFIG$DataInter, "signal.cgi.matrix.bed"),xlab="Distance to CGI Center (bp)")
plot.group.methy.profile(-4000, 4000,'SCLC',file.path(CONFIG$DataInter, "signal.cgi.matrix.bed"),xlab="Distance to CGI Center (bp)")
dev.off()
saveImage("methylation.level.profile.tss.pdf",width = 8,height = 2)
par(mfrow=c(1, 4))
plot.group.methy.profile(-4000, 4000,'LUAD',file.path(CONFIG$DataInter, "signal.tss.matrix.bed"),xlab="Distance to TSS Center (bp)")
plot.group.methy.profile(-4000, 4000,'LUSC',file.path(CONFIG$DataInter, "signal.tss.matrix.bed"),xlab="Distance to TSS Center (bp)")
plot.group.methy.profile(-4000, 4000,'LCC',file.path(CONFIG$DataInter, "signal.tss.matrix.bed"),xlab="Distance to TSS Center (bp)")
plot.group.methy.profile(-4000, 4000,'SCLC',file.path(CONFIG$DataInter, "signal.tss.matrix.bed"),xlab="Distance to TSS Center (bp)")
dev.off()

#==============================================================================================================
# Part 2
#
# The *.matrix.bed are generate by `methytools region-methy`
#==============================================================================================================
plot.region.methy<-function(signal.file){
  ratio <- loadData(file.path(CONFIG$DataInter, signal.file),header = TRUE)
  data<-data.frame(
    group=SAMPLE$sample2group(colnames(ratio)),
    ratio=unlist(ratio)
  )
  ggplot(data=data,aes(x=group,y=ratio,fill=group))+
    scale_fill_manual(values=COLOR_MAP_GROUP) +
    geom_violin( colour = 'NA',trim = TRUE)+
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="#ffffff",size=0.2)+
    theme_classic()+
    stat_compare_means( comparisons = list(c('CTL','LUAD'),c('CTL','LUSC'),c('CTL','SCLC'),c('CTL','LCC')),
                        label = 'p.signif', method = "t.test")+
    labs(x='',y='Mean Methylation Level')+
    theme(legend.position="none",
          axis.title.x = element_text(size=0),
          axis.title.y = element_text(size=14),
          axis.text = element_text(size = 14,colour="black"))+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}

saveImage('methylation.level.cgi.n50_p50.mean.pdf', width = 2.5,height = 2.5)
plot.region.methy('signal.cgi_n50_p50.matrix.bed')
dev.off()
saveImage('methylation.level.cgi.p3900_p4000.mean.pdf', width = 2.5,height = 2.5)
plot.region.methy('signal.cgi_p3900_p4000.matrix.bed')
dev.off()
saveImage('methylation.level.tss.n50_p50.mean.pdf', width = 2.5,height = 2.5)
plot.region.methy('signal.tss_n50_p50.matrix.bed')
dev.off()
saveImage('methylation.level.tss.p3900_p4000.mean.pdf', width = 2.5,height = 2.5)
plot.region.methy('signal.tss_p3900_p4000.matrix.bed')
dev.off()


#----------------------------------------------------------------------------------------------------------------------
# Figure S1: cluster of all lung cancer samples
#----------------------------------------------------------------------------------------------------------------------
methy10x<-readRDS(file.path(CONFIG$DataRaw, 'merge.d10.one.bed.rds'))
set.seed(100)
idx<-sample(1:nrow(methy10x),200000)
m<-methy10x[idx, 4:ncol(methy10x)]
dist_mm<-dist(t(m))
hclust_avg <- hclust(dist_mm,method='ward.D2')
dend <- as.dendrogram(hclust_avg)
labels_colors(dend) <- SAMPLE$table$Color[order.dendrogram(dend)]
plot(dend)

saveImage("methylation.level.cluster.pdf",width = 8,height = 8)
phyl<-as.phylo(hclust_avg)
plot(phyl, type = "fan",tip.color=SAMPLE$table$Color,label.offset=15)
tiplabels(pch=21, col="black", adj=0, bg=SAMPLE$table$Color, cex=2)
dev.off()
saveImage("methylation.level.cluster.legend.pdf",width = 5,height = 3)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(COLOR_MAP_GROUP),pch=16, pt.cex=3, cex=1.5, bty='n',col = COLOR_MAP_GROUP,horiz=TRUE)
dev.off()





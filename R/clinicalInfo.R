source('R/base.R')

library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
#----------------------------------------------------------------------------------------------------------------------
# Figure 1C: The methylation level density distribution of CTL vs. AIS, CTL vs. MIA and CTL vs. IAC
#----------------------------------------------------------------------------------------------------------------------
data<-read_excel(file.path(CONFIG$DataRaw, 'SupplementaryData.xlsx'),sheet = 'qc')
data$MCALL_MeanRatioCG_3X


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
        # axis.line = element_line(linewidth = 1),none
        axis.title.x = element_text(size=0),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size = 14,colour="black"))+
  guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1C: The methylation level density distribution of CTL vs. AIS, CTL vs. MIA and CTL vs. IAC
#----------------------------------------------------------------------------------------------------------------------

ratio <- loadData(file.path(CONFIG$DataRaw, 'merge.group.3x.all.bed'),header = TRUE)
ratio<-removeNegativeOne(ratio)



d1 <- density(ratio$CTL)
d2 <- density(ratio$LUAD)
d3 <- density(ratio$LUSQ)
d4 <- density(ratio$LCC)
d5 <- density(ratio$SCLC)
dens <- list(a=d1,b=d2,c=d3,d=d4, e=d5)
plot(NA, xlim=range(sapply(dens, "[", "x")),
     ylim=range(sapply(dens, "[", "y")),xlab="",ylab="Density",cex.lab=1.5, cex.axis=1.5)
mapply(lines, dens, col=COLOR_MAP_GROUP)
polygon(d1, col=NA,border = COLOR_MAP_GROUP[match('CTL', names(COLOR_MAP_GROUP))])
polygon(d2, col=NA, border=COLOR_MAP_GROUP[match('LUAD', names(COLOR_MAP_GROUP))])
polygon(d3, col=NA, border=COLOR_MAP_GROUP[match('LUSC', names(COLOR_MAP_GROUP))])
polygon(d4, col=NA, border=COLOR_MAP_GROUP[match('LCC', names(COLOR_MAP_GROUP))])
polygon(d5, col=NA, border=COLOR_MAP_GROUP[match('SCLC', names(COLOR_MAP_GROUP))])
#


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




signal.file<-file.path(CONFIG$DataInter, "signal.tss.matrix.bed")


x<-seq(UP, DWON)
signalMatrixGroup<-loadData(signal.file)
mmData<-signalMatrixGroup[,-1]
mm<-mmData[(10001+UP):(10001+DWON),]
# mm<-t(scale(t(mm)))
mm<-as.data.frame(apply(mm,2,smooth2(51)))
mm<-as.data.frame(mm)
par(mar = c(5,5,1,1))
if (draw.legend){
  plot(NA, xlim=c(UP, DWON), xlab=xlab,ylab=ylab,ylim=c(min(mm),max(mm)), bty='n')
  # box()
}else{
  plot(NA, xlim=c(UP, DWON),cex.lab=1.5, cex.axis=1.5,ylim=c(min(mm),max(mm)), xaxt='n',yaxt='n',xlab = "",ylab = "")
}
lines(x,mm[["CTL"]],col=g2color('CTL',alpha = 0.7))
lines(x,mm[[group]],col=g2color(group,alpha = 0.7))






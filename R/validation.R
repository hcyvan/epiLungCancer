source('R/base.R')
source('R/base.vector.R')
library(easyepi)

# Public data ------------------------------------------------------------------
## GSE186458 -------------------------------------------------------------------
### CpG ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/w4/LUAD.sample.mvm'))
data<-mvm$getByCpG(174592)
plotWindow2(data,window = 4)

### Sample ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/w4/LUAD.sample.mvm'))
data<-mvm$getBySample('Lung-Alveolar-Epithelial-Z000000T1')
plotWindow2(data,window=4)
### SPCR -----------------------------------------------------------------------
plotSPCR<-function(spcr,hide.names=TRUE){
  mvh<-MVH(spcr)
  m<-mvh$matrix()
  m<-m[rowSums(m!=0)!=0,]
  m<-t(m)
  # m <- m[order(rowSums(m), decreasing = TRUE), ]
  row_ha = rowAnnotation(bar = anno_barplot(rowSums(m)))
  m1<-log(m+0.0001)
  Heatmap(m1,
          name='Count',
          col=colorRamp2(c(0,max(m1)), c("#fffeee", "#c82423")),
          right_annotation = row_ha,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          show_row_names = !hide.names,
          show_column_names = !hide.names
  )
}

saveImage("spcr.GSE186458.luad.0.pdf",width = 20,height = 30)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LUAD_1.0_0.4.spcr'),hide.names = FALSE)
dev.off()
saveImage("spcr.GSE186458.luad.pdf",width = 6,height = 10)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LUAD_1.0_0.4.spcr'))
dev.off()
saveImage("spcr.GSE186458.lcc.0.pdf",width = 20,height = 30)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LCC_1.0_0.4.spcr'),hide.names = FALSE)
dev.off()
saveImage("spcr.GSE186458.lcc.pdf",width = 6,height = 10)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LCC_1.0_0.4.spcr'))
dev.off()
saveImage("spcr.GSE186458.lusc.0.pdf",width = 20,height = 30)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LUSC_1.0_0.4.spcr'),hide.names = FALSE)
dev.off()
saveImage("spcr.GSE186458.lusc.pdf",width = 6,height = 10)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/LUSC_1.0_0.4.spcr'))
dev.off()
saveImage("spcr.GSE186458.sclc.0.pdf",width = 20,height = 30)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/SCLC_1.0_0.4.spcr'),hide.names = FALSE)
dev.off()
saveImage("spcr.GSE186458.sclc.pdf",width = 6,height = 10)
plotSPCR(file.path(CONFIG$DataInter,'vector/GSE186458/w4/spcr/SCLC_1.0_0.4.spcr'))
dev.off()
## GSE79279 -------------------------------------------------------------------
# select more SMVs ------------------------------------------------------------
groups<-names(COLOR_MAP_GROUP)[-1]
files<-lapply(groups, function(x){
  list(
    a=file.path(CONFIG$DataInter,'vector/w4',paste0(x,'_1.0_0.2.smvc')),
    b=file.path(CONFIG$DataInter,'vector/w4',paste0(x,'.smvc'))
  )
})
names(files)<-groups
files[[3]]$a<-file.path(CONFIG$DataInter,'vector/w4','LCC_1.0_0.3.smvc')
smvcs<-lapply(files, function(x){
  SMVC(x$a,x$b)
})
saveRDS(smvcs,file.path(CONFIG$DataInter,'vector/w4','smvcs.forGSE7929.rds'))
smvcs<-readRDS(file.path(CONFIG$DataInter,'vector/w4','smvcs.forGSE7929.rds'))

### CpG ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE79279/w4/LUAD.sample.mvm'))
data<-mvm$getByCpG(14627311)
plotWindow2(data,window = 4)
### Sample ---------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE79279/w4/LUAD.sample.mvm'))
data<-mvm$getBySample('NC-P-1')
plotWindow2(data,window=4)
### SPCR -----------------------------------------------------------------------


#### LUSC ----------------------------------------------------------------------
##### plasma -------------------------------------------------------------------
mvh<-MVH(file.path(CONFIG$DataInter,'vector/GSE79279/w4/spcr/LUSC.spcr'))
m<-mvh$matrix()
a<-colSums(m)
a<-a[sapply(names(a), function(x){grepl('LC-P',x)})]
a<-sort(a,decreasing = FALSE)
typeLUSC<-c('LC-P-2','LC-P-4','LC-P-5','LC-P-7','LC-P-11','LC-P-21')
data<-data.frame(sample=names(a), hit=a, color='#666666')
data$color[match(typeLUSC, data$sample)]<-'#c82423'
data$sample<-factor(data$sample, levels = data$sample)
saveImage("spcr.GSE79279.GSE79277.plasma.lusc.pdf",width = 2,height = 3)
ggplot(data, aes(x=sample , y=hit)) +
  geom_segment( aes(xend=sample, yend=0)) +
  geom_point( size=2, color=data$color) +
  labs(x = NULL, y = NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_flip()+
  theme(
    axis.text.y = element_text(hjust = 1),
    axis.ticks.y = element_line(),
    axis.line.y.right = element_line()
  ) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right")
dev.off()

##### tissue -------------------------------------------------------------------
mvh<-MVH(file.path(CONFIG$DataInter,'vector/GSE79279/w4/spcr/LUSC.spcr'))
m<-mvh$matrix()
a<-colSums(m)
a<-a[sapply(names(a), function(x){grepl('LC-T',x)})]
a<-sort(a,decreasing = FALSE)
typeLUSC<-c('LC-T-2','LC-T-4','LC-T-5')
data<-data.frame(sample=names(a), hit=a, color='#666666')
data$color[match(typeLUSC, data$sample)]<-'#c82423'
data$sample<-factor(data$sample, levels = data$sample)
saveImage("spcr.GSE79279.GSE79277.tissue.lusc.pdf",width = 2,height = 0.9)
ggplot(data, aes(x=sample , y=hit)) +
  geom_segment( aes(xend=sample, yend=0)) +
  geom_point( size=2, color=data$color) +
  labs(x = NULL, y = NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_flip()+
  theme(
    axis.text.y = element_text(hjust = 1),
    axis.ticks.y = element_line(),
    axis.line.y.right = element_line()
  ) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right")
dev.off()
#### SCLC ----------------------------------------------------------------------
mvh<-MVH(file.path(CONFIG$DataInter,'vector/GSE79279/w4/spcr/SCLC.spcr'))
m<-mvh$matrix()
a<-colSums(m)
a<-a[sapply(names(a), function(x){grepl('N37',x)})]
a<-sort(a,decreasing = FALSE)
data<-data.frame(sample=names(a), hit=a, color='#c82423')
data$sample<-factor(data$sample, levels = data$sample)

saveImage("spcr.GSE79279.GSE79215.tissue.sclc.pdf",width = 2,height = 2)
ggplot(data, aes(x=sample , y=hit)) +
  geom_segment( aes(xend=sample, yend=0)) +
  geom_point( size=2, color=data$color) +
  labs(x = NULL, y = NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_flip()+
  theme(
    axis.text.y = element_text(hjust = 1),
    axis.ticks.y = element_line(),
    axis.line.y.right = element_line()
  ) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right")
dev.off()




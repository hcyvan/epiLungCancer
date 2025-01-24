source('R/base.R')
source('R/base.vector.R')
library(easyepi)


#########################################

uxm.file<-file.path(CONFIG$DataInter,'vector/w4.2/lungWithGSE186458.leukocytes/uxm/all.cfDNA.uxm.bed')
uxm<-read.csv(uxm.file,sep='\t',check.names = FALSE)

m<-uxm[,6:ncol(uxm)]

Heatmap(m,
        name='uxm',
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE
)

#########################################

uxm.file<-file.path(CONFIG$DataInter,'vector/w4.2/lungWithGSE186458.leukocytes/uxm/all.mock.uxm.bed')
uxm<-read.csv(uxm.file,sep='\t',check.names = FALSE)

m<-uxm[,6:ncol(uxm)]

Heatmap(m,
        name='uxm',
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE
)



#######################################################################################

all_data<-read_excel(file.path(CONFIG$DataRaw,"AllData.xlsx"),sheet = "jiashan")
uxm.file<-file.path(CONFIG$DataInter,'vector/w4.2/lungWithGSE186458.leukocytes/uxm/all.uxm.bed')
uxm<-read.csv(uxm.file,sep='\t',check.names = FALSE)
table.ori<-filter(SAMPLE$table.ori, !is.na(SampleNameTissue))
m<-uxm[, table.ori$WGBS_data]
colnames(m)<-table.ori$SampleNameTissue
m[m==-1]<-0


groups<-SAMPLE$sample2group(colnames(m))
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=groups),
  col = list(Stage =COLOR_MAP_GROUP),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
data.group<-c(rep('LUAD',500), rep('LUSC',500),rep('LCC',500),rep('SCLC',500))
row_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=factor(data.group,levels = FACTOR_LEVEL_GROUP)),
  col = list(Stage =COLOR_MAP_GROUP),
  show_annotation_name =FALSE,
  which = 'row'
)
# saveImage("uxm.top500.pdf",width = 8,height = 7)
Heatmap(m,
        name='uxm',
        col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
        top_annotation = column_annotation,
        right_annotation  = row_annotation,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE
)
# dev.off()

#############
cfDNA<-read_excel(file.path(CONFIG$DataRaw,"SupplementaryData.xlsx"),sheet = "cfDNA")
cfDNA<-cfDNA[,c('SampleName','Group','Usage')]
cfDNA<-cfDNA%>%filter(!is.na(SampleName))
samples<-cfDNA$SampleName
# m<-m[,names(sort(colSums(m!=0)))]
COLOR_MAP<-c('#2878b5','#c82423','#ffb15f')
names(COLOR_MAP)<-c('Healthy','LUAD','LUSC')


m<-uxm[1:500,6:ncol(uxm)]
m[m==-1]<-0
m<-m[,samples[samples%in%colnames(m)]]
column_annotation <-HeatmapAnnotation(
  df=data.frame(aa=cfDNA$Group),
  zz = anno_barplot(colSums(m!=0)),
  col = list(aa =COLOR_MAP),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage("uxm.cfDNA.LUAD.top500.pdf",width = 5,height = 2.5)
Heatmap(m,
        name='Count',
        top_annotation = column_annotation,
        col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE
)
dev.off()


m<-uxm[501:1000,6:ncol(uxm)]
m[m==-1]<-0
m<-m[,samples[samples%in%colnames(m!=0)]]
column_annotation <-HeatmapAnnotation(
  df=data.frame(aa=cfDNA$Group),
  zz = anno_barplot(colSums(m!=0)),
  col = list(aa =COLOR_MAP),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage("uxm.cfDNA.LUSC.top500.pdf",width = 5,height = 2.5)
Heatmap(m,
        name='Count',
        top_annotation = column_annotation,
        col=colorRamp2(c(0,max(m)), c("#fffeee", "#c82423")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE
)
dev.off()
################################################################################
Healthy.t





# Public data ------------------------------------------------------------------
## GSE186458 -------------------------------------------------------------------
### CpG ------------------------------------------------------------------------
mvm<-MVM(file.path(CONFIG$DataInter,'vector/GSE186458/w4/LUAD.sample.mvm'))
data<-mvm$getByCpG(174592)
plotWindow2(data,window = 4)

a<-rownames(data)
tissues<-sapply(a, function(x){
  xx<-strsplit(x, '-')[[1]]
  xx<-xx[1:(length(xx)-1)]
  paste(xx, collapse ='-')
})%>%unique()
length(tissues)
dim(data)
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

#### aaaa ----------------------------------------------------------------------
all_data<-read_excel(file.path(CONFIG$DataRaw,"AllData.xlsx"),sheet = "jiashan")
plasma<-filter(all_data, SampleType=='plasma')
mvh<-MVH(file.path(CONFIG$DataInter,'vector/cfDNA/w4/spcr/lungWithGSE186458.leukocytes/LUSC_0.95_0.2.spcr'))
m<-mvh$matrix()
# m<-m[rowSums(m!=0)!=0,]
# m <- m[order(rowSums(m), decreasing = TRUE), ]
# m1<-m1[,sample(1:ncol(m1),ncol(m1)/2)]
# m1<-m1[,1:ncol(m1)/2]
# m0<-m[1:2000,]

m<-m[,plasma$SampleId]
m<-m[,order(colSums(m))]

m0<-log(m+0.0001)
# m0<-m
m0<-m0[sample(1:nrow(m0), 500),]
Heatmap(m0,
        name='Count',
        # col=colorRamp2(c(0,max(m0)), c("#fffeee", "#c82423")),
        cluster_rows = FALSE,
        cluster_columns = FALSE
        # show_row_names = !hide.names,
        # show_column_names = !hide.names
)

#### aaaa ----------------------------------------------------------------------



m[,plasma$SampleId]%>%dim


match(plasma$SampleId, colnames(m))

source('R/base.R')

library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)

data<-read_excel(file.path(CONFIG$DataRaw, 'SupplementaryData.xlsx'),sheet = 'qc')
data$MCALL_MeanRatioCG_3X


samplesMatch<-data[match(SAMPLE$table$SampleName,data$SampleName),]

data<-data.frame(
  group=samples$Group,
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
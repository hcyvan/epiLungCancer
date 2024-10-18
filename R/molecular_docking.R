source('R/base.R')
library(easyepi)



md<-read.csv(file.path(CONFIG$DataInter,'md/md.csv'))
md$Disease<-factor(md$Disease, levels = c('NSCLC','SCLC','Depression','SCZ','AD','PD'))
md<-arrange(md,md$Disease)

data1<-data.frame(x=md$Drug, y=md$AverageBindingEnergy, disease=md$Disease,color='#c82423')
data2<-data.frame(x=md$Drug, y=md$DominantBindingEnergy,disease=md$Disease ,color='#2878b5')
# data3<-data.frame(x=md$Drug, y=md$Average_CDOCKER_ENERGY,disease=md$Disease ,color='green3')
# data4<-data.frame(x=md$Drug, y=md$Dominant_CDOCKER_ENERGY,disease=md$Disease ,color='orange')
data<-rbind(data1,data2)
data$x<-factor(data$x, levels =as.character(md$Drug))
saveImage("md.result.pdf",width = 2.5,height = 4)
ggplot(data, aes(x=x , y=y)) +
  geom_segment(aes(xend=x, yend=0)) +
  geom_point( size=1.6, color=data$color) +
  labs(x = NULL, y = NULL)+
  geom_hline(yintercept = -6,color = "gray")+
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
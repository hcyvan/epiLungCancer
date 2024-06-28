library(tidyverse)

#reading
dir<-''
raw <- read.csv(paste0(dir,"raw.csv"),row.names = 1)
#           GSM1218807 GSM1218808 GSM1218809 GSM1218810 GSM1218811 GSM1218812 GSM1218813
#cg00000292 0.44405010 0.43441214 0.51796360 0.59439646 0.10782704 0.49950950 0.34169978
#cg00002426 0.40140440 0.47621528 0.36323052 0.51230966 0.18715256 0.39787500 0.58386209
cg <- read.csv(paste0(dir,"cg.csv"))#Annotated information from the GPL (GEO database)
# ID CHR   MAPINFO chrom     start       end
# cg00000292  16  28890100   chr16  28890099  28890101
# cg00002426   3  57743543    chr3  57743542  57743544
# cg00003994   7  15725862    chr7  15725861  15725863
# cg00006414   7 148822837    chr7 148822836 148822838
# cg00007981  11  93862594   chr11  93862593  93862595
# cg00008493  14  93813777   chr14  93813776  93813778


#match ID
df<-cbind(rownames(raw),raw)
colnames(df)[1]<-'ID'
a = cg[match(df$ID, cg$ID),]
a[,7]<-paste0(a[,4],':',a[,3],'-',a[,6])
df[,1]<-a[,7]
colnames(df)[1]<-'site'
write.csv(a[,c(4:6)],paste0(dir,"hg19.csv"),row.names = F)

#hg19-hg38 transformation by hgLiftOver (Convert to bed format)

#match hg38 information
b <- read.csv(paste0(dir,"hg38.csv"))
x<-b[match(df$site, b$site),]
y<-cbind(x[,1:3],df[,-1])
y<-y %>%
  filter(chrom != '')

write.csv(y,paste0(dir,"raw_hg38.csv"),row.names = F)








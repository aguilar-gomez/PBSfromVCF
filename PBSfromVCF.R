library(ggplot2)
library(latex2exp)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
#Tile
custom_title<-args[1]
#Population name 0
pbs0<-args[2]
#Population name 1
pbs1<-args[3]
#Population name 2
pbs2<-args[4]
#For wiritng outlier file
outquantile<-args[5]
#Color 1
color1<-args[6]
#File concating Fst 
persitefst <- read.delim(args[7], header=TRUE)
#name for outfile
name<-args[8]

p0 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST))
p1 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST.1))
p2 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST.2))
persitefst$pbs0 = (p0 + p1 - p2) / 2
persitefst$pbs1 = (p0 + p2 - p1) / 2
persitefst$pbs2 = (p1 + p2 - p0) / 2

#Get range of PBS
gezero<-persitefst[persitefst$pbs0>0,]
pbs_fil<-gezero[gezero$pbs0<10,]

write.table(pbs_fil[c("CHROM","POS","pbs0","pbs1","pbs2")],
            paste0("pbs_",name,".tsv"),
            quote = F,sep = "\t",row.names = F)

#Only plot 25% of the data
toplot<-quantile(pbs_fil$pbs0,.75)
pbs_fil2plot<-pbs_fil[pbs_fil$pbs0>toplot,]

#Find outliers
cutoff<-quantile(pbs_fil2plot$pbs0,as.numeric(outquantile))
hits<-pbs_fil2plot[pbs_fil2plot$pbs0>=cutoff,]
hits$CHROM<-paste0("chr",hits$CHROM)
write.table(hits[c("CHROM","POS","POS","pbs0")],
            paste0(pbs0,"_pbs_outliers_top_persite",outquantile,".bed"),
            quote = F,sep = "\t",row.names = F)

pbs_fil2plot$SNPs<-1:nrow(pbs_fil2plot)
pbs_fil2plot$CHROM<-paste0("chr",pbs_fil2plot$CHROM)

axis_set <- pbs_fil2plot %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))

nolabel<-ggplot(pbs_fil2plot, aes(x=SNPs,y=pbs0,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "gray"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  xlab("Chromosome")+ylab(TeX(paste0("$PBS_{",pbs0,"}$")))  


ggsave(paste0(pbs0,"_pbs0_persite.png"),width=10,height=5)


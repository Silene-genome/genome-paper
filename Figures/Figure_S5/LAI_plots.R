library(ggplot2)

t <- read.table("S.latifolia_v4.0.genome_v2.fna.out.LAI",sep = "\t",header=T)
head(t)

chrOrder<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chrX","chrY")
t$Chr<-factor(t$Chr,levels=chrOrder)

chrCols<-c(rep(c("#08519c","#4292c6"),5),"#08519c","#41ab5d","#cc4c02")

laiPlot<-ggplot(t, aes(x=Chr,y=LAI,fill=Chr)) + 
  geom_boxplot() + 
  geom_hline(yintercept=20.38, col='red',size=1.5) + 
  scale_fill_manual(values=chrCols) + 
  ggtitle("LTR assambly index [LAI] in S. latifolia",subtitle="whole genome LAI = 20.38 [horizontal red line]") + 
  ylab("LAI") + xlab("Chromosome") + 
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),legend.position="None")

# point plot
png("Slat_LAI_points.png",res=300,units="in",width=15, height=8)
ggplot(t, aes(x=Chr,y=LAI,fill=Chr)) + 
  geom_point(aes(colour=factor(Chr),
             fill=factor(Chr)),
             size=3,
             position = position_dodge2(w = 1),
             alpha=0.5) + 
  geom_hline(yintercept=20.38, col='red',size=1.5) + 
  scale_color_manual(values=chrCols) + 
  #ggtitle("LTR assambly index [LAI] in S. latifolia",subtitle="whole genome LAI = 20.38 [horizontal red line]") + 
  ylab("LAI score") + xlab("Chromosome") + 
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),legend.position="None")
dev.off()

# chrSpec
png("Slat_LAI_pointsChrSpec.png",res=300,units="in",width=10, height=15)
ggplot(t, aes(x=From_Mbp,y=LAI,fill=Chr)) + 
  geom_point(aes(colour=factor(Chr),
                 fill=factor(Chr)),
             alpha=0.5) + 
  geom_hline(yintercept=20.38, col='red',size=1) + 
  scale_color_manual(values=chrCols) + 
  #ggtitle("LTR assambly index [LAI] in S. latifolia",subtitle="whole genome LAI = 20.38 [horizontal red line]") + 
  ylab("LAI score") + xlab("Chromosome") + 
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.position="None")+ 
  facet_grid(rows=vars(Chr))
dev.off()


png("Slat_LAI.png",res=300,units="in",width=8, height=7)
laiPlot
dev.off()

library(scales)
t$From_Mbp <- t$From/1000000

png("Slat_LAI_chrSpecific.png",res=300,units="in",width=8, height=15)
ggplot(t, aes(x=From_Mbp,y=LAI,fill=Chr)) + 
  geom_bar(stat='identity') + 
  geom_hline(yintercept=20.38, col='red',size=1.5) + 
  scale_fill_manual(values=chrCols) + 
  ggtitle("LTR assambly index [LAI] in S. latifolia",subtitle="whole genome LAI = 20.38 [horizontal red line]") + 
  ylab("LAI") + xlab("Sliding window position FROM [Mbp]") + 
  scale_x_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.position="None")+ 
  facet_grid(rows=vars(Chr))

dev.off()
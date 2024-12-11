# SupplFig March 31 2024

library(ggplot2)
#ageCols_v2<-c('#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850')
ageCols_v2<-rev(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'))

#ageCols_v2<-rev(c('#ffffcc','#ffffcc','#fed976','#fed976','#fd8d3c','#fd8d3c','#e31a1c','#e31a1c','#e31a1c'))
chrOrder<-c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chrX',
            'chrY')

t_ac <- read.table("chrY_MYA_chDistr_ageCat_chrYMod.tab",sep="\t")
colnames(t_ac)<-c("ch", "fam","sec","start", "end", "edta_ltrSim", "ltr_sim", "ltr_avg", "ltr_l", "ltr_r", "K80", "K2P_loc", "K2P_glob", "mya","bps","age_cat")
t_ac$ch<- factor(t_ac$ch,levels = chrOrder)

t_ac_famSub <- subset(t_ac,t_ac$fam=='copia_Ale'|t_ac$fam=='copia_Angela'|t_ac$fam=='copia_SIRE'|t_ac$fam=='gypsy_Athila'|t_ac$fam=='gypsy_CRM'|t_ac$fam=='gypsy_Ogre'|t_ac$fam=='gypsy_Retand'|t_ac$fam=='gypsy_Tat'|t_ac$fam=='gypsy_Tekay')


#png("Slat_v4_TE_ageDistribution_chr3_April16.png",res=300,units = 'in',width = 10,height = 10)

ggplot(subset(t_ac_famSub,t_ac_famSub$ch=="chr3"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  ggtitle("LTR retrotransposon age distribution in chr3 of S. latifolia genome",
          subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TE proportion [TE family+age category bps per chromosome bps]") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=13)) +
  facet_grid(rows = vars(fam))
dev.off()



png("Slat_v4_TE_ageDistribution_chrX_April16.png",res=300,units = 'in',width = 10,height = 10)

ggplot(subset(t_ac_famSub,t_ac_famSub$ch=="chrX"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  ggtitle("LTR retrotransposon age distribution in chrX of S. latifolia genome",
          subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TE proportion [TE family+age category bps per chromosome bps]") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=13)) +
  facet_grid(rows = vars(fam))
dev.off()


png("Slat_v4_TE_ageDistribution_chrY_April16.png",res=300,units = 'in',width = 10,height = 10)

ggplot(subset(t_ac_famSub,t_ac_famSub$ch=="chrY"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  ggtitle("LTR retrotransposon age distribution in chrY of S. latifolia genome",
          subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TE proportion [TE family+age category bps per chromosome bps]") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=13)) +
  facet_grid(rows = vars(fam))
dev.off()


head(t_ac_famSub)
# Define the background colors for each facet
background_colors <- c(rep("yellow", 3), rep("blue", 6))

# Define the background colors for row and column facets
row_background_colors <- c("yellow", "yellow", "yellow")
column_background_colors <- c("lightblue", "lightblue", "lightblue", "lightblue","lightblue", "lightblue")

pLeg <-ggplot(subset(t_ac_famSub, t_ac_famSub$ch == "chrY"), aes(x = sec, y = bps, fill = age_cat)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values = ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels = c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                               '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  ggtitle("LTR retrotransposon age distribution in chrY of S. latifolia genome",
          subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TE proportion [TE family+age category bps per chromosome bps]") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 13),
        strip.background = element_rect(fill = row_background_colors),
        strip.background.x = element_rect(fill = column_background_colors)) +  # Set background colors for row and column facets
  facet_grid(rows = vars(fam))

t_ac_famSub$fam_copy <- t_ac_famSub$fam


library(tidyr)
# Split the original column by "_"
t_ac_famSub_sep <- separate(t_ac_famSub, fam, into = c("supFam", "fFam"), sep = "_")
head(t_ac_famSub_sep)
fFamOrder <-c('Ale','Angela','SIRE','Athila','CRM','Ogre','Retand','Tat','Tekay')

t_ac_famSub_sep$fFam <- factor(t_ac_famSub_sep$fFam,levels = fFamOrder)
t_ac_copia <- subset(t_ac_famSub_sep, t_ac_famSub_sep$supFam == 'copia')
t_ac_gypsy <- subset(t_ac_famSub_sep, t_ac_famSub_sep$supFam == 'gypsy')


library(ggplot2)
library(tidyverse)
library(grid)

p <- ggplot(subset(t_ac_famSub_sep,t_ac_famSub_sep$ch=="chrY"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) +
  
  #ggtitle("LTR retrotransposon age distribution in chrY of S. latifolia genome",
  #        subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TEs proportion [bps per Mbp]") +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=11),
        strip.text = element_text(color = "white",size=13),
        legend.position = 'None') +
  facet_grid(rows = vars(fFam))

#p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))

fills <- c(rep("#66bd63", 3), rep("#4393c3", 6))
# fills <- c("red","green","blue","yellow")
 k <- 1

 for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
 }
 

grid.draw(g)

library(ggplot2)
p1 <- ggplot(subset(t_ac_famSub_sep,t_ac_famSub_sep$ch=="chrX"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) +
  
  #ggtitle("LTR retrotransposon age distribution in chrY of S. latifolia genome",
  #        subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TEs proportion [bps per Mbp]") +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=11),
        strip.text = element_text(color = "white",size=13),
        legend.position = 'None') +
  facet_grid(rows = vars(fFam))

#p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl)
g1 <- ggplot_gtable(ggplot_build(p1))
stripr <- which(grepl('strip-r', g1$layout$name))

fills <- c(rep("#66bd63", 3), rep("#4393c3", 6))
# fills <- c("red","green","blue","yellow")
k <- 1

for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


grid.draw(g1)


p2 <- ggplot(subset(t_ac_famSub_sep,t_ac_famSub_sep$ch=="chr3"),aes(x=sec,y=bps,fill=age_cat)) +
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(name = 'TE age category [MYA]:',
                    values=ageCols_v2,
                    breaks = c('A','B','C','D','E','F','G','H','I'),
                    labels=c('0.0-0.5', '0.5-1.0', '1.0-1.5','1.5-2.0','2.0-2.5',
                             '2.5-3.0','3.0-3.5','3.5-4.0','> 4.0')) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) +
  #ggtitle("LTR retrotransposon age distribution in chrY of S. latifolia genome",
  #        subtitle = "Intact TEs only [DANTE and EDTA combination]") +
  xlab("Chromosome section [Mbp]") +
  ylab("TEs proportion [bps per Mbp]") +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=11),
        strip.text = element_text(color = "white",size=13),
        legend.position = 'None') +
  facet_grid(rows = vars(fFam))

#p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl)
g2 <- ggplot_gtable(ggplot_build(p2))
stripr <- which(grepl('strip-r', g2$layout$name))

fills <- c(rep("#66bd63", 3), rep("#4393c3", 6))
# fills <- c("red","green","blue","yellow")
k <- 1

for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g2)

png("chY_TEage.png",res=300,units = 'in',width = 8,height = 8)
grid.draw(g)
dev.off()

png("chX_TEage.png",res=300,units = 'in',width = 8,height = 8)
grid.draw(g1)
dev.off()


png("ch3_TEage.png",res=300,units = 'in',width = 8,height = 8)
grid.draw(g2)
dev.off()

library(patchwork)
(p+p1+p2)

library(cowplot)
legend <- cowplot::get_legend(pLeg)

grid.newpage()
grid.draw(legend)
ggsave("legend_plot.png", plot = legend, width = 2, height = 5, units = "in")


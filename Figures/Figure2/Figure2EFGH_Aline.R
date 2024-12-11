#!/usr/bin/env Rscript


library(data.table)
library(ggplot2)


DirectoryPath <- "/home/muyle/Documents/RECHERCHE/DNA_methylation_Silene/Silene_latifolia_genome_paper"
setwd(DirectoryPath)


#---------------------------
# Y/X ration (Figure 2.d.1)
#---------------------------
XYpairs <- fread("XYpairs_gene_expr_sRNA.txt", h=T)

pdf("Figure_2.d.1_Y_over_X_hist.pdf", width=5, height=5)
ggplot(aes(x = Y_over_X_gene_expression), data=XYpairs[tissue=="flowers"]) +
  theme_bw(base_size=20) +
  geom_histogram(aes(x = Y_over_X_gene_expression), bins=3000, color = "black", size = 0.2) +
  geom_boxplot(aes(x = Y_over_X_gene_expression, y=-1), fun="mean", show.legend = FALSE, col="red3", outlier.shape = NA) +
  geom_vline(xintercept = 1, col = "black", linetype = 3) +
  coord_cartesian(xlim=c(0, 3)) +
  labs(x="Y over X gene expression in males", y="Number of X/Y gene pairs") #+
#  facet_grid(rows = vars(tissue), cols = vars(strata), scales = "free_y")
dev.off()  

#-----------------------------
# 24nt sRNA RPM (Figure 2.d.3)
#-----------------------------
XYpairs_Alleles <- fread("XYpairs_Alleles_sRNA.txt", h=T)

pdf("Figure_2.d.3_24_sRNA_mapping.pdf", width=5, height=5)
ggplot(aes(y = sRNA24_exp, x = Allele, fill=Allele), data=XYpairs_Alleles[tissue=="flowers"]) +
  theme_bw(base_size=20) +
  geom_boxplot(aes(y = sRNA24_exp, x = Allele, fill=Allele), show.legend = FALSE, outlier.shape = NA) +
  stat_summary(aes(y = sRNA24_exp, x = Allele), fun="mean", show.legend = FALSE, col="red3") +
  coord_cartesian(ylim=c(0, 0.2)) +
  labs(y="24nt sRNA mapping (RPM)", x="Allele") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) #+
  #  facet_grid(rows = vars(tissue), scales = "free_y")
#  facet_grid(rows = vars(tissue), cols = vars(strata), scales = "free_y")
dev.off()  


#---------------------------------------------------------------
# CHH DNA methylation windows for X/Y gene pairs (Figure 2.d.4)
#---------------------------------------------------------------
gene_SummaryWindows <- fread("gene_SummaryWindows.txt", h=T)

gene_SummaryWindows$Allele <- factor(gene_SummaryWindows$Allele, levels = c('XX in female', 'X in male', 'Y in male'))

pdf("Figure_2.d.4_methylation_XYpairs.pdf", width=7, height=5)
ggplot(data=gene_SummaryWindows[context=='CHH'], aes(x = window, y = Average_Methylation, col=Allele)) + 
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_point(aes(x = window, y = Average_Methylation, col=Allele), size = 1) +  
  geom_line(aes(x = window, y = Average_Methylation, col=Allele), size=0.5) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high, col=Allele, fill=Allele), linetype=1, size=0.2, alpha=0.1, show.legend = FALSE) +
  labs(x=paste("position", sep=''), y="X/Y genes CHH methylation") + 
  guides(col=guide_legend(title="Allele:")) +
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000')) #+
#  facet_grid(rows = vars(context), scales = "free_y")
dev.off()


#----------------------------------------------------
# TEs +/- 4000bp around XY gene pairs (Figure 2.d.2)
#----------------------------------------------------
XY_genes_repeats <- fread('XY_genes_repeats.txt', h=T)

pdf("Figure_2.d.2_TE_length_around_XY_pairs.pdf", width=5, height=5)
ggplot(data=XY_genes_repeats, aes(x =chr, y = total_length_TEs, fill = chr)) + 
  theme_bw(base_size=20) +
  geom_boxplot(aes(y = total_length_TEs, fill = chr), outlier.shape = NA, show.legend = FALSE) +
  stat_summary(aes(y = total_length_TEs, x=chr), fun="mean", show.legend = FALSE, col="red3") +
  labs(x="Allele", y="Repeat length around X/Y genes", fill="Allele") +
  coord_cartesian(ylim=c(0, 35000))# +
#  facet_grid(rows = vars(strata), scales = "free_y")
dev.off()

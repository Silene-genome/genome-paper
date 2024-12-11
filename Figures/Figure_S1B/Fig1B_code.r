# Datapath in line 25 needs to be changed based on the file location.

# Add fig lab function - from Sophie/web ####
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {
  
  text(x = line2user(xl, 2), y = line2user(yl, 3), 
       lab, xpd = NA, font = 2, cex = 1.5, adj = c(0, 1))
  
}

# Colour functions ####
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}

# for DE plots
myred <- t_col("red", 80)
myblue <- t_col("blue", 80)


## Input data ####

###  Marey plot ####
datapath <- "~/Desktop/Silene_Code/Figures/Fig1B/rec rate"
m_data <- read.table(file.path(datapath, "male_loess.txt"), header=T)
m_data$recombinationRate[m_data$recombinationRate<0] <- 0
f_data <- read.table(file.path(datapath, "female_loess.txt"), header=T)
f_data$recombinationRate[f_data$recombinationRate<0] <- 0


### F/M coverage
datapath <- "~/Desktop/Silene_Code/Figures/Fig1B/FM ratio-sex linked genes-circos"
c_data <- read.table(file.path(datapath, "v4_depth_male_female.txt"), header=T)


### differential expression ####

datapath <- "~/Desktop/Silene_Code/Figures/Fig1B/deg"
deg5_data <- read.csv(file.path(datapath, "table_gene_DE_S5_F-S5_M.with_annotation.csv"), header=T)
deg8_data <- read.csv(file.path(datapath, "table_gene_DE_S8_F-S8_M.with_annotation.csv"), header=T)


outpath <- "~/Desktop/Silene_Code/Figures/Fig1B"



## MAIN FIGURE

par(mfrow=c(3,1), xpd=T, oma=c(4,0,0,0), mar=c(0,6,2,2))

x_start=290

# cM ####
plot(c(x_start, 350), c(0,120), xlab="Mb", type="n", ylab="", cex.lab=1, las=1, bty="n", xaxt="n")
Axis(side=1, labels=FALSE)
lines(m_data$phys[m_data$map=="chrX" & m_data$phys>289000000]/1000000, m_data$gen[m_data$map=="chrX" & m_data$phys>289000000], lw=3, col="blue")
lines(f_data$phys[f_data$map=="chrX" & f_data$phys>289000000]/1000000, f_data$gen[f_data$map=="chrX" & f_data$phys>289000000], lw=3, col="red")
title(ylab="cM", line=4, cex.lab=1.5)
segments(x0=320224904/1000000, x1=320224904/1000000, y0=0, y1=115, lty=2, lwd=2)
legend("topleft", inset=c(0.02,0.05),
       legend=c("male","female"),
       col=c('blue',"red"), ncol=1, cex=1, bty="n", lwd=2)
text(x=321, y=115, "PAR", adj=0, font=2, cex=1.5)
addfiglab("B")


# F/M coverage ####
plot(c(x_start, 350), c(0,3), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n", xaxt="n")
Axis(side=1, labels=FALSE)
lines(c_data$start[c_data$chr=="chrX" & c_data$start>290000000]/1000000, c_data$ratioF_M[c_data$chr=="chrX" & c_data$start>290000000], lwd=3, col="black")
title(ylab=paste("F/M \n coverage", sep=""), line=3, cex.lab=1.5)
segments(x0=x_start, x1=346, y0=1, y1=1, lty=2)
segments(x0=x_start, x1=346, y0=2, y1=2, lty=2)
segments(x0=320224904/1000000, x1=320224904/1000000, y0=0, y1=3, lty=2, lwd=2)


# differential expression ####
plot(c(x_start, 350), c(-10,10), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n", xaxt="n")
Axis(side=1)
points(deg5_data$beg[deg5_data$chr=="chr12" & deg5_data$log2FC>0 & deg5_data$beg>290000000]/1000000, deg5_data$log2FC[deg5_data$chr=="chr12" & deg5_data$log2FC>0 & deg5_data$beg>290000000], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr12" & deg5_data$log2FC < 0 & deg5_data$beg>290000000]/1000000, deg5_data$log2FC[deg5_data$chr=="chr12" & deg5_data$log2FC < 0 & deg5_data$beg>290000000], lw=2, pch=20, col=myblue)
title(ylab=paste("DE genes \n stage 5", sep=""), line=3, cex.lab=1.5)
abline(v=320224904/1000000, lty=2, lwd=2, xpd=F)
segments(x0=x_start, x1=348, y0=0, y1=0, lty=2)
legend("topleft", inset=c(0.0,0.05),
       legend=c("male-biased","female-biased"),
       col=c(myblue,myred), ncol=1, cex=1, bty="n", pch=20, xpd=T)

mtext(text="Mb", side=1, line=2.5, outer=TRUE)
quartz.save(type = 'pdf', file = paste(outpath, '/Fig1B.pdf', sep=""))
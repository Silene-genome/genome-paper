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

myred <- t_col("red", 80)
myblue <- t_col("blue", 80)


datapath <- "~/Desktop/Silene_Code/Figures/FigS18"
deg5_data <- read.csv(file.path(datapath, "table_gene_DE_S5_F-S5_M.with_annotation.csv"), header=T)
deg8_data <- read.csv(file.path(datapath, "table_gene_DE_S8_F-S8_M.with_annotation.csv"), header=T)



# Panel A
par(mfrow=c(6,2),mar=c(4,5.5,2,1), xpd=T, oma=c(1,0,0,0))

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr1"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr1" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr1" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr1" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr1" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr1", line=4, cex.lab=1)
addfiglab("A")

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr2"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr2" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr2" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr2" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr2" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr2", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr3"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr3" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr3" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr3" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr3" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr3", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr4"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr4" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr4" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr4" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr4" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr4", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr5"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr5" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr5" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr5" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr5" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr5", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr6"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr6" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr6" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr6" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr6" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr6", line=4, cex.lab=1)
ree
plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr7"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr7" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr7" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr7" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr7" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr7", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr8"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr8" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr8" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr8" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr8" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr8", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr9"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr9" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr9" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr9" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr9" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr9", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr10"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr10" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr10" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr10" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr10" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr10", line=4, cex.lab=1)

plot(c(0, max(deg5_data$beg[deg5_data$chr=="chr11"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr11" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr11" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr11" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr11" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr11", line=4, cex.lab=1)


plot(c(0, 340), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg5_data$beg[deg5_data$chr=="chr12" & deg5_data$log2FC>0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr12" & deg5_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg5_data$beg[deg5_data$chr=="chr12" & deg5_data$log2FC < 0]/1000000, deg5_data$log2FC[deg5_data$chr=="chr12" & deg5_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr12", line=4, cex.lab=1)
abline(v=320224904/1000000, lty=2, lwd=2, xpd=F)
legend("top", inset=c(0,-0.1),
       legend=c("male-biased","female-biased"),
       col=c('blue',"red"), ncol=2, cex=.9, bty="n", lwd=2, xpd=T)


quartz.save(type = 'pdf', file = paste(datapath, '/FigS18A.pdf', sep=""))



# Panel B
par(mfrow=c(6,2),mar=c(4,5.5,2,1), xpd=T, oma=c(1,0,0,0))

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr1"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr1" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr1" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr1" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr1" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr1", line=4, cex.lab=1)

addfiglab("B")

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr2"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr2" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr2" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr2" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr2" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr2", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr3"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr3" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr3" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr3" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr3" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr3", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr4"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr4" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr4" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr4" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr4" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr4", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr5"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr5" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr5" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr5" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr5" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr5", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr6"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr6" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr6" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr6" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr6" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr6", line=4, cex.lab=1)
ree
plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr7"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr7" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr7" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr7" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr7" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr7", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr8"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr8" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr8" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr8" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr8" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr8", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr9"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr9" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr9" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr9" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr9" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr9", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr10"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr10" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr10" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr10" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr10" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr10", line=4, cex.lab=1)

plot(c(0, max(deg8_data$beg[deg8_data$chr=="chr11"])/1000000), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr11" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr11" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr11" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr11" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr11", line=4, cex.lab=1)


plot(c(0, 340), c(-18,14), xlab="", type="n", ylab="", cex.lab=1, las=1, bty="n")
points(deg8_data$beg[deg8_data$chr=="chr12" & deg8_data$log2FC>0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr12" & deg8_data$log2FC>0], lw=2, pch=20, col=myred)
points(deg8_data$beg[deg8_data$chr=="chr12" & deg8_data$log2FC < 0]/1000000, deg8_data$log2FC[deg8_data$chr=="chr12" & deg8_data$log2FC < 0], lw=2, pch=20, col=myblue)
title(ylab="chr12", line=4, cex.lab=1)
abline(v=320224904/1000000, lty=2, lwd=2, xpd=F)
legend("top", inset=c(0,-0.1),
       legend=c("male-biased","female-biased"),
       col=c('blue',"red"), ncol=2, cex=.9, bty="n", lwd=2, xpd=T)


quartz.save(type = 'pdf', file = paste(datapath, '/FigS18B.pdf', sep=""))
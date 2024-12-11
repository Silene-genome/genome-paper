# Silene analysis autosomes, X and Y
# 2023-09-12
# Sophie
#set directory ####
setwd("~/Documents/PROJECTS/Silene consortium/Data/Figure S14")

##load packages #### 
library(qqman)
library(dplyr)
library(zoo)
library(ggplot2)
library(gridExtra)                                         # Install scales package
library(scales)
library(ggridges)
library(gridExtra)

# Add fig lab function####
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








#FST within-species FRAPOL####----------------------------------
# read September 4 2023 data####

getwd()
###read data from txt####
FRA1POL3_fst <- read.delim("female/fst/FRA1POL3_fst.txt")
str(FRA1POL3_fst)
table(FRA1POL3_fst$chromosome)

chrY.FRA1POL3.haploid_fst <- read.delim("male/Fst/chrY.FRA1POL3.haploid_fst.txt")
str(chrY.FRA1POL3.haploid_fst)
table(chrY.FRA1POL3.haploid_fst$chromosome)

FRA1POL3_fst <- rbind(FRA1POL3_fst, chrY.FRA1POL3.haploid_fst)
###mid window length & mid position####
FRA1POL3_fst$window_length <- FRA1POL3_fst$window_pos_2- FRA1POL3_fst$window_pos_1
FRA1POL3_fst$mid_pos <- FRA1POL3_fst$window_pos_1 + FRA1POL3_fst$window_length/2


###remove Fst < 0, Fst =NA, window length < 100 and adjust chromosome and SNP var ####
table(FRA1POL3_fst$avg_hudson_fst <0)
table(FRA1POL3_fst$window_length < 100)
table(is.na(FRA1POL3_fst$avg_hudson_fst))
table(FRA1POL3_fst$avg_hudson_fst <0 | FRA1POL3_fst$window_length < 100 | is.na(FRA1POL3_fst$avg_hudson_fst)==T)



FRA1POL3_fst_m <- FRA1POL3_fst %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(avg_hudson_fst <0), NA)) %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(window_length < 100), NA)) %>% 
  mutate(SNP=seq(1,dim(FRA1POL3_fst)[1], by=1), .after=avg_hudson_fst) %>% 
  na.omit() 

dim(FRA1POL3_fst)[1] - dim(FRA1POL3_fst_m)[1]
str(FRA1POL3_fst_m)

#define chromosome types
FRA1POL3_fst_m$chr.type <- "Autosomes"
FRA1POL3_fst_m$chr.type[FRA1POL3_fst_m$chromosome=="chrX"] <- "X-chr."
FRA1POL3_fst_m$chr.type[FRA1POL3_fst_m$chromosome=="chrY"] <- "Y-chr."
table(FRA1POL3_fst_m$chr.type)
FRA1POL3_fst_m$chr.type <- factor(FRA1POL3_fst_m$chr.type)

###distance between SNP + stats ####

FRA1POL3_fst_m$dist <- NA

for(n in 2:dim(FRA1POL3_fst_m)[1]) 
{
  FRA1POL3_fst_m$dist[n]  <- FRA1POL3_fst_m$window_pos_1[n]-FRA1POL3_fst_m$window_pos_1[n-1]
}

sum(FRA1POL3_fst_m$dist <0, na.rm=T)
FRA1POL3_fst_m$dist[FRA1POL3_fst_m$dist<0] <- NA

summary(FRA1POL3_fst_m$dist)
stats.all_FRAPOL <- data.frame(comparison="within_FRAPOL", chr.type=levels(FRA1POL3_fst_m$chr.type), 
                               no.RADtags.tot = NA, 
                               no.SNP.tot = NA,
                               no.SNP.per.RADtag.med=NA, 
                               no.SNP.per.RADtag.q1=NA, 
                               no.SNP.per.RADtag.q3=NA,
                               window_length.med=NA,
                               window_length.min=NA,
                               window_length.max=NA,
                               dist.med=NA,
                               dist.q1=NA,
                               dist.q3=NA
)

stats.all_FRAPOL$no.RADtags.tot <- table(FRA1POL3_fst_m$chr.type)
stats.all_FRAPOL$no.SNP.tot <- tapply(FRA1POL3_fst_m$no_snps, INDEX=FRA1POL3_fst_m$chr.type, 
                                      FUN=function(x) sum(x, na.rm=T))
stats.all_FRAPOL$no.SNP.per.RADtag.med <- round(tapply(FRA1POL3_fst_m$no_snps, INDEX=FRA1POL3_fst_m$chr.type, 
                                                       FUN=function(x) median(x, na.rm=T)), 1)

stats.all_FRAPOL$no.SNP.per.RADtag.q1 <- round(tapply(FRA1POL3_fst_m$no_snps, INDEX=FRA1POL3_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.25, na.rm=T)), 1)

stats.all_FRAPOL$no.SNP.per.RADtag.q3 <- round(tapply(FRA1POL3_fst_m$no_snps, INDEX=FRA1POL3_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.75, na.rm=T)), 1)

stats.all_FRAPOL$window_length.med <- round(tapply(FRA1POL3_fst_m$window_length, INDEX=FRA1POL3_fst_m$chr.type, 
                                                   FUN=function(x) median(x, na.rm=T)), 1)

stats.all_FRAPOL$window_length.min <- round(tapply(FRA1POL3_fst_m$window_length, INDEX=FRA1POL3_fst_m$chr.type, 
                                                   FUN=function(x) min(x, na.rm=T)),1)

stats.all_FRAPOL$window_length.max <- round(tapply(FRA1POL3_fst_m$window_length, INDEX=FRA1POL3_fst_m$chr.type, 
                                                   FUN=function(x) max(x, na.rm=T)),1)

stats.all_FRAPOL$dist.med <- round(tapply(FRA1POL3_fst_m$dist, INDEX=FRA1POL3_fst_m$chr.type, 
                                          FUN=function(x) median(x, na.rm=T)), 0)
stats.all_FRAPOL$dist.q1 <- round(tapply(FRA1POL3_fst_m$dist, INDEX=FRA1POL3_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.25, na.rm=T)),0)
stats.all_FRAPOL$dist.q3 <- round(tapply(FRA1POL3_fst_m$dist, INDEX=FRA1POL3_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.75, na.rm=T)),0)


write.table(stats.all_FRAPOL, "Stats output/stats.all_FRAPOL_230912.csv", sep=",")



###rolling median and rolling mean####
#fix chromosome numbering
str(FRA1POL3_fst_m)
FRA1POL3_fst_m$chromosome.orig <- FRA1POL3_fst_m$chromosome
FRA1POL3_fst_m$chromosome <- gsub("[chr]", "",FRA1POL3_fst_m$chromosome.orig)
table(FRA1POL3_fst_m$chromosome)
FRA1POL3_fst_m$chromosome[FRA1POL3_fst_m$chromosome=="X"] <- 12
FRA1POL3_fst_m$chromosome[FRA1POL3_fst_m$chromosome=="Y"] <- 13

a <- 12
for(n in 1:a){
  FRA1POL3_fst_m$movmed.51[FRA1POL3_fst_m$chromosome==n] <- 
    rollmedian( FRA1POL3_fst_m$avg_hudson_fst[FRA1POL3_fst_m$chromosome==n], 51, fill=NA) 
  FRA1POL3_fst_m$movmean.50[FRA1POL3_fst_m$chromosome==n] <- 
    rollmean( FRA1POL3_fst_m$avg_hudson_fst[FRA1POL3_fst_m$chromosome==n], 50, fill=NA) 
}

str(FRA1POL3_fst_m)

###set start positions for Manhattan plot ####


a<- 13
start <- data.frame(chromosome=c(1:a), max=1:a, plotstart=1:a, min=1:a, length=NA)

for(n in 1:a){
  start$max[n] <- max(FRA1POL3_fst_m$window_pos_1[FRA1POL3_fst_m$chromosome==n], na.rm=T)
}

a
for(n in 1:a){
  start$min[n] <- min(FRA1POL3_fst_m$window_pos_1[FRA1POL3_fst_m$chromosome==n], na.rm=T)
}

start$plotstart[1]<- 0
#from Xiadong 23-09-07
#chr1	200709468
#chr2	199635652
#chr3	153359320
#chr4	112764170
#chr5	91147008
#chr6	141279031
#chr7	128518071
#chr8	113857337
#chr9	226323612
#chr10	165488795
#chr11	209616475
#chrX	346484273
#chrY	497814031

start$plotstart[1]<- 0
start$length <- c(200709468,
                  199635652,
                  153359320,
                112764170,
                  91147008,
                  141279031,
                  128518071,
                  113857337,
                  226323612,
                  165488795,
                209616475,
                  346484273,
                497814031)

for(n in 2:a){
  start$plotstart[n] <- sum(start$length[start$chromosome<n], na.rm=T)
}

start$tick <- NA
for (n in 1:a){
  start$tick[n] <- start$plotstart[n]+start$length[n]/2 
}

###add windowed calculations####
#non-overlapping windows
window.size <- 5000000
a<- 12
n <- a
for (n in c(1:a))
{
  no <- floor(start$length[n]/window.size)
  current <- data.frame(chr = as.numeric(paste(n)), window=1:no, fst=NA, mid=NA, start=NA, snp=NA, scaff=NA, seq.len=NA)
  shift <- (start$length[n]-window.size*no)/2
  current$start <- seq(from=shift, to = shift + window.size*no-1, by=window.size)
  current$mid <- current$start + window.size/2
  chr <- FRA1POL3_fst_m[FRA1POL3_fst_m$chromosome==n,] 
  m <- 1
  for (m in c(1:no)){
    win <-     chr[chr$window_pos_1 >= current$start[m] 
                   & chr$window_pos_2 < current$start[m]+window.size,]
    current$fst[m] <- mean(win$avg_hudson_fst)
    current$scaff[m] <- dim(win)[1]
    current$snp[m] <- sum(win$no_snps)
    current$seq.len[m] <- sum(win$window_length)
    assign(paste("chr",n, sep="_"), current)
  }
  
}

data.wind_FRAPOL <- rbind(chr_1, chr_2, chr_3, chr_4, chr_5, chr_6, chr_7, chr_8, chr_9, chr_10, chr_11, chr_12)


data.wind_FRAPOL$chr.type[data.wind_FRAPOL$chr < 12] <- "Autosomes"
data.wind_FRAPOL$chr.type[data.wind_FRAPOL$chr==12] <- "X-chr."

table(data.wind_FRAPOL$chr.type)
data.wind_FRAPOL$chr.type <- factor(data.wind_FRAPOL$chr.type)



sum(data.wind_FRAPOL$scaff < 4, na.rm = T)
data.wind_FRAPOL <- data.wind_FRAPOL[data.wind_FRAPOL$scaff > 4,]

window.summary_FRAPOL <- data.frame(chr.type=levels(data.wind_FRAPOL$chr.type), window.size=window.size,
                                    no.RADtags.per.wind.med =NA, 
                                    no.RADtags.q1 =NA, 
                                    no.RADtags.q3 =NA,
                                    no.SNP.per.wind.med=NA,
                                    no.SNP.per.wind.q1=NA, 
                                    no.SNP.per.wind.q3=NA,
                                    seq.length.per.wind.med=NA,
                                    seq.length.per.wind.q1=NA,
                                    seq.length.per.wind.q3=NA)

window.summary_FRAPOL$no.RADtags.per.wind.med <- round(tapply(data.wind_FRAPOL$scaff, data.wind_FRAPOL$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_FRAPOL$no.RADtags.q1 <- round(tapply(data.wind_FRAPOL$scaff, data.wind_FRAPOL$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_FRAPOL$no.RADtags.q3 <- round(tapply(data.wind_FRAPOL$scaff, data.wind_FRAPOL$chr.type, function(x) quantile(x, 0.75, na.rm=T)),1)

window.summary_FRAPOL$no.SNP.per.wind.med <- round(tapply(data.wind_FRAPOL$snp , data.wind_FRAPOL$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_FRAPOL$no.SNP.per.wind.q1 <- round(tapply(data.wind_FRAPOL$snp , data.wind_FRAPOL$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_FRAPOL$no.SNP.per.wind.q3 <- round(tapply(data.wind_FRAPOL$snp , data.wind_FRAPOL$chr.type, function(x) quantile(x, 0.75, na.rm=T)), 1)


window.summary_FRAPOL$seq.length.per.wind.med <- round(tapply(data.wind_FRAPOL$seq.len, data.wind_FRAPOL$chr.type, function(x) median(x, na.rm=T)),0)
window.summary_FRAPOL$seq.length.per.wind.q1 <- round(tapply(data.wind_FRAPOL$seq.len, data.wind_FRAPOL$chr.type, function(x) quantile(x, 0.25,na.rm=T)),0)
window.summary_FRAPOL$seq.length.per.wind.q3 <- round(tapply(data.wind_FRAPOL$seq.len, data.wind_FRAPOL$chr.type, function(x) quantile(x,0.75, na.rm=T)),0)

window.summary_FRAPOL
write.csv(window.summary_FRAPOL, "Stats output/betw_FRAPOL_window_stats_230912.csv")
tail(data.wind_FRAPOL)

###add start & end of  chromosome ####
#adjust no column as needed
str(FRA1POL3_fst_m)

#from Xiadong 23-09-07
#chr1	200709468
#chr2	199635652
#chr3	153359320
#chr4	112764170
#chr5	91147008
#chr6	141279031
#chr7	128518071
#chr8	113857337
#chr9	226323612
#chr10	165488795
#chr11	209616475
#chrX	346484273
#chrY	497814031

matrix(nrow=26, ncol=15, NA) -> add.end

add.end[,c(3,10)] <- matrix(c(10,	165488795,
                              11,	209616475,
                              12,	346484273,
                              13,497814031,
                              1,	200709468,
                              2,	199635652,
                              3,	153359320,
                              4,	112764170,
                              5,	91147008,
                              6,	141279031,
                              7,	128518071,
                              8,	113857337,
                              9,	226323612,
                              10,	0,
                              11,	0,
                              12,	0,
                              13,0,
                              1,	0,
                              2,	0,
                              3,	0,
                              4,	0,
                              5,	0,
                              6,	0,
                              7,	0,
                              8,	0,
                              9,	0), byrow=T, ncol=2, nrow=26)

add.end <- data.frame(add.end)

add.end[,7]<-15
colnames(add.end) <- colnames(FRA1POL3_fst_m)
str(add.end)

rbind(FRA1POL3_fst_m, add.end)-> FRA1POL3_fst_m.x
str(FRA1POL3_fst_m.x)

str(FRA1POL3_fst_m)







#FST between-species BELFRA####----------------------------------
# read September 4 2023 data####

###read data from txt####
BEL1FRA1_fst <- read.delim("female/fst/BEL1FRA1_fst.txt")
str(BEL1FRA1_fst)
table(BEL1FRA1_fst$chromosome)

chrY.BEL1FRA1.haploid_fst <- read.delim("male/Fst/chrY.BEL1FRA1.haploid_fst.txt")
str(chrY.BEL1FRA1.haploid_fst)
table(chrY.BEL1FRA1.haploid_fst$chromosome)

BEL1FRA1_fst <- rbind(BEL1FRA1_fst, chrY.BEL1FRA1.haploid_fst)
###mid window length & mid position####
BEL1FRA1_fst$window_length <- BEL1FRA1_fst$window_pos_2- BEL1FRA1_fst$window_pos_1
BEL1FRA1_fst$mid_pos <- BEL1FRA1_fst$window_pos_1 + BEL1FRA1_fst$window_length/2


###remove Fst < 0, Fst =NA, window length < 100 and adjust chromosome and SNP var ####
table(BEL1FRA1_fst$avg_hudson_fst <0)
table(BEL1FRA1_fst$window_length < 100)
table(is.na(BEL1FRA1_fst$avg_hudson_fst))
table(BEL1FRA1_fst$avg_hudson_fst <0 | BEL1FRA1_fst$window_length < 100 | is.na(BEL1FRA1_fst$avg_hudson_fst)==T)



BEL1FRA1_fst_m <- BEL1FRA1_fst %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(avg_hudson_fst <0), NA)) %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(window_length < 100), NA)) %>% 
  mutate(SNP=seq(1,dim(BEL1FRA1_fst)[1], by=1), .after=avg_hudson_fst) %>% 
  na.omit() 

dim(BEL1FRA1_fst)[1] - dim(BEL1FRA1_fst_m)[1]
str(BEL1FRA1_fst_m)

#define chromosome types
BEL1FRA1_fst_m$chr.type <- "Autosomes"
BEL1FRA1_fst_m$chr.type[BEL1FRA1_fst_m$chromosome=="chrX"] <- "X-chr."
BEL1FRA1_fst_m$chr.type[BEL1FRA1_fst_m$chromosome=="chrY"] <- "Y-chr."
table(BEL1FRA1_fst_m$chr.type)
BEL1FRA1_fst_m$chr.type <- factor(BEL1FRA1_fst_m$chr.type)

###distance between SNP + stats ####

BEL1FRA1_fst_m$dist <- NA

for(n in 2:dim(BEL1FRA1_fst_m)[1]) 
{
  BEL1FRA1_fst_m$dist[n]  <- BEL1FRA1_fst_m$window_pos_1[n]-BEL1FRA1_fst_m$window_pos_1[n-1]
}

sum(BEL1FRA1_fst_m$dist <0, na.rm=T)
BEL1FRA1_fst_m$dist[BEL1FRA1_fst_m$dist<0] <- NA

summary(BEL1FRA1_fst_m$dist)
stats.all_BELFRA <- data.frame(comparison="between_BELFRA", chr.type=levels(BEL1FRA1_fst_m$chr.type), 
                               no.RADtags.tot = NA, 
                               no.SNP.tot = NA,
                               no.SNP.per.RADtag.med=NA, 
                               no.SNP.per.RADtag.q1=NA, 
                               no.SNP.per.RADtag.q3=NA,
                               window_length.med=NA,
                               window_length.min=NA,
                               window_length.max=NA,
                               dist.med=NA,
                               dist.q1=NA,
                               dist.q3=NA
)

stats.all_BELFRA$no.RADtags.tot <- table(BEL1FRA1_fst_m$chr.type)
stats.all_BELFRA$no.SNP.tot <- tapply(BEL1FRA1_fst_m$no_snps, INDEX=BEL1FRA1_fst_m$chr.type, 
                                      FUN=function(x) sum(x, na.rm=T))
stats.all_BELFRA$no.SNP.per.RADtag.med <- round(tapply(BEL1FRA1_fst_m$no_snps, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                       FUN=function(x) median(x, na.rm=T)), 1)

stats.all_BELFRA$no.SNP.per.RADtag.q1 <- round(tapply(BEL1FRA1_fst_m$no_snps, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.25, na.rm=T)), 1)

stats.all_BELFRA$no.SNP.per.RADtag.q3 <- round(tapply(BEL1FRA1_fst_m$no_snps, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.75, na.rm=T)), 1)

stats.all_BELFRA$window_length.med <- round(tapply(BEL1FRA1_fst_m$window_length, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                   FUN=function(x) median(x, na.rm=T)), 1)

stats.all_BELFRA$window_length.min <- round(tapply(BEL1FRA1_fst_m$window_length, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                   FUN=function(x) min(x, na.rm=T)),1)

stats.all_BELFRA$window_length.max <- round(tapply(BEL1FRA1_fst_m$window_length, INDEX=BEL1FRA1_fst_m$chr.type, 
                                                   FUN=function(x) max(x, na.rm=T)),1)

stats.all_BELFRA$dist.med <- round(tapply(BEL1FRA1_fst_m$dist, INDEX=BEL1FRA1_fst_m$chr.type, 
                                          FUN=function(x) median(x, na.rm=T)), 0)
stats.all_BELFRA$dist.q1 <- round(tapply(BEL1FRA1_fst_m$dist, INDEX=BEL1FRA1_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.25, na.rm=T)),0)
stats.all_BELFRA$dist.q3 <- round(tapply(BEL1FRA1_fst_m$dist, INDEX=BEL1FRA1_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.75, na.rm=T)),0)


write.table(stats.all_BELFRA, "Stats output/stats.all_BELFRA_230912.csv", sep=",")



###rolling median and rolling mean####
#fix chromosome numbering
str(BEL1FRA1_fst_m)
BEL1FRA1_fst_m$chromosome.orig <- BEL1FRA1_fst_m$chromosome
BEL1FRA1_fst_m$chromosome <- gsub("[chr]", "",BEL1FRA1_fst_m$chromosome.orig)
table(BEL1FRA1_fst_m$chromosome)
BEL1FRA1_fst_m$chromosome[BEL1FRA1_fst_m$chromosome=="X"] <- 12
BEL1FRA1_fst_m$chromosome[BEL1FRA1_fst_m$chromosome=="Y"] <- 13

a <- 12
for(n in 1:a){
  BEL1FRA1_fst_m$movmed.51[BEL1FRA1_fst_m$chromosome==n] <- 
    rollmedian( BEL1FRA1_fst_m$avg_hudson_fst[BEL1FRA1_fst_m$chromosome==n], 51, fill=NA) 
  BEL1FRA1_fst_m$movmean.50[BEL1FRA1_fst_m$chromosome==n] <- 
    rollmean( BEL1FRA1_fst_m$avg_hudson_fst[BEL1FRA1_fst_m$chromosome==n], 50, fill=NA) 
}

str(BEL1FRA1_fst_m)

###set start positions for Manhattan plot ####


a<- 13
start <- data.frame(chromosome=c(1:a), plotstart=1:a, length=NA)



start$plotstart[1]<- 0
#from Xiadong 23-09-07
#chr1	200709468
#chr2	199635652
#chr3	153359320
#chr4	112764170
#chr5	91147008
#chr6	141279031
#chr7	128518071
#chr8	113857337
#chr9	226323612
#chr10	165488795
#chr11	209616475
#chrX	346484273
#chrY	497814031

start$plotstart[1]<- 0
start$length <- c(200709468,
                  199635652,
                  153359320,
                  112764170,
                  91147008,
                  141279031,
                  128518071,
                  113857337,
                  226323612,
                  165488795,
                  209616475,
                  346484273,
                  497814031)

for(n in 2:a){
  start$plotstart[n] <- sum(start$length[start$chromosome<n], na.rm=T)
}

start$tick <- NA
for (n in 1:a){
  start$tick[n] <- start$plotstart[n]+start$length[n]/2 
}
start
###add windowed calculations####
#non-overlapping windows
window.size <- 5000000
a<- 12
n <- a
for (n in c(1:a))
{
  no <- floor(start$length[n]/window.size)
  current <- data.frame(chr = as.numeric(paste(n)), window=1:no, fst=NA, mid=NA, start=NA, snp=NA, scaff=NA, seq.len=NA)
  shift <- (start$length[n]-window.size*no)/2
  current$start <- seq(from=shift, to = shift + window.size*no-1, by=window.size)
  current$mid <- current$start + window.size/2
  chr <- BEL1FRA1_fst_m[BEL1FRA1_fst_m$chromosome==n,] 
  m <- 1
  for (m in c(1:no)){
    win <-     chr[chr$window_pos_1 >= current$start[m] 
                   & chr$window_pos_2 < current$start[m]+window.size,]
    current$fst[m] <- mean(win$avg_hudson_fst)
    current$scaff[m] <- dim(win)[1]
    current$snp[m] <- sum(win$no_snps)
    current$seq.len[m] <- sum(win$window_length)
    assign(paste("chr",n, sep="_"), current)
  }
  
}

data.wind_BELFRA <- rbind(chr_1, chr_2, chr_3, chr_4, chr_5, chr_6, chr_7, chr_8, chr_9, chr_10, chr_11, chr_12)


data.wind_BELFRA$chr.type[data.wind_BELFRA$chr < 12] <- "Autosomes"
data.wind_BELFRA$chr.type[data.wind_BELFRA$chr==12] <- "X-chr."

table(data.wind_BELFRA$chr.type)
data.wind_BELFRA$chr.type <- factor(data.wind_BELFRA$chr.type)



sum(data.wind_BELFRA$scaff < 4, na.rm = T)
data.wind_BELFRA <- data.wind_BELFRA[data.wind_BELFRA$scaff > 4,]

window.summary_BELFRA <- data.frame(chr.type=levels(data.wind_BELFRA$chr.type), window.size=window.size,
                                    no.RADtags.per.wind.med =NA, 
                                    no.RADtags.q1 =NA, 
                                    no.RADtags.q3 =NA,
                                    no.SNP.per.wind.med=NA,
                                    no.SNP.per.wind.q1=NA, 
                                    no.SNP.per.wind.q3=NA,
                                    seq.length.per.wind.med=NA,
                                    seq.length.per.wind.q1=NA,
                                    seq.length.per.wind.q3=NA)

window.summary_BELFRA$no.RADtags.per.wind.med <- round(tapply(data.wind_BELFRA$scaff, data.wind_BELFRA$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_BELFRA$no.RADtags.q1 <- round(tapply(data.wind_BELFRA$scaff, data.wind_BELFRA$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_BELFRA$no.RADtags.q3 <- round(tapply(data.wind_BELFRA$scaff, data.wind_BELFRA$chr.type, function(x) quantile(x, 0.75, na.rm=T)),1)

window.summary_BELFRA$no.SNP.per.wind.med <- round(tapply(data.wind_BELFRA$snp , data.wind_BELFRA$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_BELFRA$no.SNP.per.wind.q1 <- round(tapply(data.wind_BELFRA$snp , data.wind_BELFRA$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_BELFRA$no.SNP.per.wind.q3 <- round(tapply(data.wind_BELFRA$snp , data.wind_BELFRA$chr.type, function(x) quantile(x, 0.75, na.rm=T)), 1)


window.summary_BELFRA$seq.length.per.wind.med <- round(tapply(data.wind_BELFRA$seq.len, data.wind_BELFRA$chr.type, function(x) median(x, na.rm=T)),0)
window.summary_BELFRA$seq.length.per.wind.q1 <- round(tapply(data.wind_BELFRA$seq.len, data.wind_BELFRA$chr.type, function(x) quantile(x, 0.25,na.rm=T)),0)
window.summary_BELFRA$seq.length.per.wind.q3 <- round(tapply(data.wind_BELFRA$seq.len, data.wind_BELFRA$chr.type, function(x) quantile(x,0.75, na.rm=T)),0)

window.summary_BELFRA
write.csv(window.summary_BELFRA, "Stats output/betw_BELFRA_window_stats_230912.csv")
tail(data.wind_BELFRA)
window.summary_FRAPOL

###add start & end of  chromosome ####
#adjust no column as needed
str(BEL1FRA1_fst_m)

#from Xiadong 23-09-07
#chr1	200709468
#chr2	199635652
#chr3	153359320
#chr4	112764170
#chr5	91147008
#chr6	141279031
#chr7	128518071
#chr8	113857337
#chr9	226323612
#chr10	165488795
#chr11	209616475
#chrX	346484273
#chrY	497814031

matrix(nrow=26, ncol=15, NA) -> add.end

add.end[,c(3,10)] <- matrix(c(10,	165488795,
                              11,	209616475,
                              12,	346484273,
                              13,497814031,
                              1,	200709468,
                              2,	199635652,
                              3,	153359320,
                              4,	112764170,
                              5,	91147008,
                              6,	141279031,
                              7,	128518071,
                              8,	113857337,
                              9,	226323612,
                              10,	0,
                              11,	0,
                              12,	0,
                              13,0,
                              1,	0,
                              2,	0,
                              3,	0,
                              4,	0,
                              5,	0,
                              6,	0,
                              7,	0,
                              8,	0,
                              9,	0), byrow=T, ncol=2, nrow=26)

add.end <- data.frame(add.end)

add.end[,7]<-15
colnames(add.end) <- colnames(BEL1FRA1_fst_m)
str(add.end)

rbind(BEL1FRA1_fst_m, add.end)-> BEL1FRA1_fst_m.x
str(BEL1FRA1_fst_m.x)

str(BEL1FRA1_fst_m)













#FST between-species BELPOL####----------------------------------
# read September 4 2023 data####


###read data from txt####
BEL1POL3_fst <- read.delim("female/fst/BEL1POL3_fst.txt")
str(BEL1POL3_fst)
table(BEL1POL3_fst$chromosome)

chrY.BEL1POL3.haploid_fst <- read.delim("male/Fst/chrY.BEL1POL3.haploid_fst.txt")
str(chrY.BEL1POL3.haploid_fst)
table(chrY.BEL1POL3.haploid_fst$chromosome)

BEL1POL3_fst <- rbind(BEL1POL3_fst, chrY.BEL1POL3.haploid_fst)
###mid window length & mid position####
BEL1POL3_fst$window_length <- BEL1POL3_fst$window_pos_2- BEL1POL3_fst$window_pos_1
BEL1POL3_fst$mid_pos <- BEL1POL3_fst$window_pos_1 + BEL1POL3_fst$window_length/2


###remove Fst < 0, Fst =NA, window length < 100 and adjust chromosome and SNP var ####
table(BEL1POL3_fst$avg_hudson_fst <0)
table(BEL1POL3_fst$window_length < 100)
table(is.na(BEL1POL3_fst$avg_hudson_fst))
table(BEL1POL3_fst$avg_hudson_fst <0 | BEL1POL3_fst$window_length < 100 | is.na(BEL1POL3_fst$avg_hudson_fst)==T)



BEL1POL3_fst_m <- BEL1POL3_fst %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(avg_hudson_fst <0), NA)) %>% 
  mutate(avg_hudson_fst = replace(avg_hudson_fst, which(window_length < 100), NA)) %>% 
  mutate(SNP=seq(1,dim(BEL1POL3_fst)[1], by=1), .after=avg_hudson_fst) %>% 
  na.omit() 

dim(BEL1POL3_fst)[1] - dim(BEL1POL3_fst_m)[1]
str(BEL1POL3_fst_m)

#define chromosome types
BEL1POL3_fst_m$chr.type <- "Autosomes"
BEL1POL3_fst_m$chr.type[BEL1POL3_fst_m$chromosome=="chrX"] <- "X-chr."
BEL1POL3_fst_m$chr.type[BEL1POL3_fst_m$chromosome=="chrY"] <- "Y-chr."
table(BEL1POL3_fst_m$chr.type)
BEL1POL3_fst_m$chr.type <- factor(BEL1POL3_fst_m$chr.type)

###distance between SNP + stats ####

BEL1POL3_fst_m$dist <- NA

for(n in 2:dim(BEL1POL3_fst_m)[1]) 
{
  BEL1POL3_fst_m$dist[n]  <- BEL1POL3_fst_m$window_pos_1[n]-BEL1POL3_fst_m$window_pos_1[n-1]
}

sum(BEL1POL3_fst_m$dist <0, na.rm=T)
BEL1POL3_fst_m$dist[BEL1POL3_fst_m$dist<0] <- NA

summary(BEL1POL3_fst_m$dist)
stats.all_BELPOL <- data.frame(comparison="between_BELPOL", chr.type=levels(BEL1POL3_fst_m$chr.type), 
                               no.RADtags.tot = NA, 
                               no.SNP.tot = NA,
                               no.SNP.per.RADtag.med=NA, 
                               no.SNP.per.RADtag.q1=NA, 
                               no.SNP.per.RADtag.q3=NA,
                               window_length.med=NA,
                               window_length.min=NA,
                               window_length.max=NA,
                               dist.med=NA,
                               dist.q1=NA,
                               dist.q3=NA
)

stats.all_BELPOL$no.RADtags.tot <- table(BEL1POL3_fst_m$chr.type)
stats.all_BELPOL$no.SNP.tot <- tapply(BEL1POL3_fst_m$no_snps, INDEX=BEL1POL3_fst_m$chr.type, 
                                      FUN=function(x) sum(x, na.rm=T))
stats.all_BELPOL$no.SNP.per.RADtag.med <- round(tapply(BEL1POL3_fst_m$no_snps, INDEX=BEL1POL3_fst_m$chr.type, 
                                                       FUN=function(x) median(x, na.rm=T)), 1)

stats.all_BELPOL$no.SNP.per.RADtag.q1 <- round(tapply(BEL1POL3_fst_m$no_snps, INDEX=BEL1POL3_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.25, na.rm=T)), 1)

stats.all_BELPOL$no.SNP.per.RADtag.q3 <- round(tapply(BEL1POL3_fst_m$no_snps, INDEX=BEL1POL3_fst_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.75, na.rm=T)), 1)

stats.all_BELPOL$window_length.med <- round(tapply(BEL1POL3_fst_m$window_length, INDEX=BEL1POL3_fst_m$chr.type, 
                                                   FUN=function(x) median(x, na.rm=T)), 1)

stats.all_BELPOL$window_length.min <- round(tapply(BEL1POL3_fst_m$window_length, INDEX=BEL1POL3_fst_m$chr.type, 
                                                   FUN=function(x) min(x, na.rm=T)),1)

stats.all_BELPOL$window_length.max <- round(tapply(BEL1POL3_fst_m$window_length, INDEX=BEL1POL3_fst_m$chr.type, 
                                                   FUN=function(x) max(x, na.rm=T)),1)

stats.all_BELPOL$dist.med <- round(tapply(BEL1POL3_fst_m$dist, INDEX=BEL1POL3_fst_m$chr.type, 
                                          FUN=function(x) median(x, na.rm=T)), 0)
stats.all_BELPOL$dist.q1 <- round(tapply(BEL1POL3_fst_m$dist, INDEX=BEL1POL3_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.25, na.rm=T)),0)
stats.all_BELPOL$dist.q3 <- round(tapply(BEL1POL3_fst_m$dist, INDEX=BEL1POL3_fst_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.75, na.rm=T)),0)


write.table(stats.all_BELPOL, "Stats output/stats.all_BELPOL_230912.csv", sep=",")



###rolling median and rolling mean####
#fix chromosome numbering
str(BEL1POL3_fst_m)
BEL1POL3_fst_m$chromosome.orig <- BEL1POL3_fst_m$chromosome
BEL1POL3_fst_m$chromosome <- gsub("[chr]", "",BEL1POL3_fst_m$chromosome.orig)
table(BEL1POL3_fst_m$chromosome)
BEL1POL3_fst_m$chromosome[BEL1POL3_fst_m$chromosome=="X"] <- 12
BEL1POL3_fst_m$chromosome[BEL1POL3_fst_m$chromosome=="Y"] <- 13

a <- 12
for(n in 1:a){
  BEL1POL3_fst_m$movmed.51[BEL1POL3_fst_m$chromosome==n] <- 
    rollmedian( BEL1POL3_fst_m$avg_hudson_fst[BEL1POL3_fst_m$chromosome==n], 51, fill=NA) 
  BEL1POL3_fst_m$movmean.50[BEL1POL3_fst_m$chromosome==n] <- 
    rollmean( BEL1POL3_fst_m$avg_hudson_fst[BEL1POL3_fst_m$chromosome==n], 50, fill=NA) 
}

str(BEL1POL3_fst_m)


###add windowed calculations####
#non-overlapping windows
window.size <- 5000000
a<- 12
n <- a
for (n in c(1:a))
{
  no <- floor(start$length[n]/window.size)
  current <- data.frame(chr = as.numeric(paste(n)), window=1:no, fst=NA, mid=NA, start=NA, snp=NA, scaff=NA, seq.len=NA)
  shift <- (start$length[n]-window.size*no)/2
  current$start <- seq(from=shift, to = shift + window.size*no-1, by=window.size)
  current$mid <- current$start + window.size/2
  chr <- BEL1POL3_fst_m[BEL1POL3_fst_m$chromosome==n,] 
  m <- 1
  for (m in c(1:no)){
    win <-     chr[chr$window_pos_1 >= current$start[m] 
                   & chr$window_pos_2 < current$start[m]+window.size,]
    current$fst[m] <- mean(win$avg_hudson_fst)
    current$scaff[m] <- dim(win)[1]
    current$snp[m] <- sum(win$no_snps)
    current$seq.len[m] <- sum(win$window_length)
    assign(paste("chr",n, sep="_"), current)
  }
  
}

data.wind_BELPOL <- rbind(chr_1, chr_2, chr_3, chr_4, chr_5, chr_6, chr_7, chr_8, chr_9, chr_10, chr_11, chr_12)


data.wind_BELPOL$chr.type[data.wind_BELPOL$chr < 12] <- "Autosomes"
data.wind_BELPOL$chr.type[data.wind_BELPOL$chr==12] <- "X-chr."

table(data.wind_BELPOL$chr.type)
data.wind_BELPOL$chr.type <- factor(data.wind_BELPOL$chr.type)



sum(data.wind_BELPOL$scaff < 4, na.rm = T)
data.wind_BELPOL <- data.wind_BELPOL[data.wind_BELPOL$scaff > 4,]

window.summary_BELPOL <- data.frame(chr.type=levels(data.wind_BELPOL$chr.type), window.size=window.size,
                                    no.RADtags.per.wind.med =NA, 
                                    no.RADtags.q1 =NA, 
                                    no.RADtags.q3 =NA,
                                    no.SNP.per.wind.med=NA,
                                    no.SNP.per.wind.q1=NA, 
                                    no.SNP.per.wind.q3=NA,
                                    seq.length.per.wind.med=NA,
                                    seq.length.per.wind.q1=NA,
                                    seq.length.per.wind.q3=NA)

window.summary_BELPOL$no.RADtags.per.wind.med <- round(tapply(data.wind_BELPOL$scaff, data.wind_BELPOL$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_BELPOL$no.RADtags.q1 <- round(tapply(data.wind_BELPOL$scaff, data.wind_BELPOL$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_BELPOL$no.RADtags.q3 <- round(tapply(data.wind_BELPOL$scaff, data.wind_BELPOL$chr.type, function(x) quantile(x, 0.75, na.rm=T)),1)

window.summary_BELPOL$no.SNP.per.wind.med <- round(tapply(data.wind_BELPOL$snp , data.wind_BELPOL$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_BELPOL$no.SNP.per.wind.q1 <- round(tapply(data.wind_BELPOL$snp , data.wind_BELPOL$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_BELPOL$no.SNP.per.wind.q3 <- round(tapply(data.wind_BELPOL$snp , data.wind_BELPOL$chr.type, function(x) quantile(x, 0.75, na.rm=T)), 1)


window.summary_BELPOL$seq.length.per.wind.med <- round(tapply(data.wind_BELPOL$seq.len, data.wind_BELPOL$chr.type, function(x) median(x, na.rm=T)),0)
window.summary_BELPOL$seq.length.per.wind.q1 <- round(tapply(data.wind_BELPOL$seq.len, data.wind_BELPOL$chr.type, function(x) quantile(x, 0.25,na.rm=T)),0)
window.summary_BELPOL$seq.length.per.wind.q3 <- round(tapply(data.wind_BELPOL$seq.len, data.wind_BELPOL$chr.type, function(x) quantile(x,0.75, na.rm=T)),0)

window.summary_BELPOL
write.csv(window.summary_BELPOL, "Stats output/betw_BELPOL_window_stats_230912.csv")
tail(data.wind_BELPOL)
window.summary_FRAPOL

###add start & end of  chromosome ####
#adjust no column as needed
str(BEL1POL3_fst_m)

#from Xiadong 23-09-07
#chr1	200709468
#chr2	199635652
#chr3	153359320
#chr4	112764170
#chr5	91147008
#chr6	141279031
#chr7	128518071
#chr8	113857337
#chr9	226323612
#chr10	165488795
#chr11	209616475
#chrX	346484273
#chrY	497814031

matrix(nrow=26, ncol=15, NA) -> add.end

add.end[,c(3,10)] <- matrix(c(10,	165488795,
                              11,	209616475,
                              12,	346484273,
                              13,497814031,
                              1,	200709468,
                              2,	199635652,
                              3,	153359320,
                              4,	112764170,
                              5,	91147008,
                              6,	141279031,
                              7,	128518071,
                              8,	113857337,
                              9,	226323612,
                              10,	0,
                              11,	0,
                              12,	0,
                              13,0,
                              1,	0,
                              2,	0,
                              3,	0,
                              4,	0,
                              5,	0,
                              6,	0,
                              7,	0,
                              8,	0,
                              9,	0), byrow=T, ncol=2, nrow=26)

add.end <- data.frame(add.end)

add.end[,7]<-15
colnames(add.end) <- colnames(BEL1POL3_fst_m)
str(add.end)

rbind(BEL1POL3_fst_m, add.end)-> BEL1POL3_fst_m.x
str(BEL1POL3_fst_m.x)

str(BEL1POL3_fst_m)















#Diversity Pi ####-------------------------------------------------------
###load and combine pi data from Xiaodong ####

###read data from txt####
FRA1POL3_pi <- read.delim("female/pi/FRA1POL3_pi.txt")
chrY.FRA1POL3.haploid_pi <- read.delim("male/pi/chrY.FRA1POL3.haploid_pi.txt")

str(FRA1POL3_pi)
str(chrY.FRA1POL3.haploid_pi)
table(chrY.FRA1POL3.haploid_pi$chromosome)

FRA1POL3_pi$chromosome[FRA1POL3_pi$chromosome=="chrX"] <- 12
chrY.FRA1POL3.haploid_pi$chromosome[chrY.FRA1POL3.haploid_pi$chromosome=="chrY"] <- 13

FRA1POL3_pi <- rbind(FRA1POL3_pi, chrY.FRA1POL3.haploid_pi)

FRA1_pi <- FRA1POL3_pi[FRA1POL3_pi$pop=="FRA1",]
POL3_pi <- FRA1POL3_pi[FRA1POL3_pi$pop=="POL3",]


str(FRA1_pi[FRA1_pi$chromosome==13,])
summary(FRA1_pi[FRA1_pi$chromosome==13,])
str(POL3_pi)

str(POL3_pi[POL3_pi$chromosome==13,])

summary(POL3_pi[POL3_pi$chromosome==13,])

#

# pi FRA1####
###mid window length & mid position####
FRA1_pi$window_length <- FRA1_pi$window_pos_2- FRA1_pi$window_pos_1
FRA1_pi$mid_pos <- FRA1_pi$window_pos_1 + FRA1_pi$window_length/2

###remove pi < 0, pi =NA, window length < 100 and adjust chromosome and SNP var ####
table(FRA1_pi$avg_pi <0)
table(FRA1_pi$window_length < 100, FRA1_pi$chromosome)
summary(FRA1_pi)
boxplot(FRA1_pi$window_length~FRA1_pi$chromosome)
table(is.na(FRA1_pi$avg_pi))
table(FRA1_pi$avg_pi <0 | FRA1_pi$window_length < 100 | is.na(FRA1_pi$avg_pi)==T)

str(FRA1_pi)

FRA1_pi_m <- FRA1_pi %>% 
  mutate(chromosome= as.numeric(gsub("chr", "", chromosome)))%>% 
  mutate(avg_pi = replace(avg_pi, which(avg_pi <0), NA)) %>% 
  mutate(avg_pi = replace(avg_pi, which(window_length < 100), NA)) %>% 
  mutate(SNP=seq(1,dim(FRA1_pi)[1], by=1), .after=avg_pi) %>% 
  na.omit() 

dim(FRA1_pi)[1] - dim(FRA1_pi_m)[1]
str(FRA1_pi_m)

#define chromosome types
FRA1_pi_m$chr.type[FRA1_pi_m$chromosome<12] <- "Autosomes"
FRA1_pi_m$chr.type[FRA1_pi_m$chromosome==12] <- "X-chr."
FRA1_pi_m$chr.type[FRA1_pi_m$chromosome==13] <- "Y-chr."

FRA1_pi_m$chr.type <- factor(FRA1_pi_m$chr.type)

###distance between SNP +stats ####

FRA1_pi_m$dist <- NA

for(n in 2:dim(FRA1_pi_m)[1]) 
{
  FRA1_pi_m$dist[n]  <- FRA1_pi_m$window_pos_1[n]-FRA1_pi_m$window_pos_1[n-1]
}

sum(FRA1_pi_m$dist <0, na.rm=T)
FRA1_pi_m$dist[FRA1_pi_m$dist<0] <- NA

summary(FRA1_pi_m$dist)
stats.all_FRA1 <- data.frame(comparison="FRA1",chr.type=levels(FRA1_pi_m$chr.type), 
                        no.RADtags.tot = NA, 
                        no.SNP.tot = NA,
                        no.SNP.per.RADtag.med=NA, 
                        no.SNP.per.RADtag.q1=NA, 
                        no.SNP.per.RADtag.q3=NA,
                        window_length.med=NA,
                        window_length.min=NA,
                        window_length.max=NA,
                        dist.med=NA,
                        dist.q1=NA,
                        dist.q3=NA
)

stats.all_FRA1$no.RADtags.tot <- table(FRA1_pi_m$chr.type)
stats.all_FRA1$no.SNP.tot <- tapply(FRA1_pi_m$no_sites, INDEX=FRA1_pi_m$chr.type, 
                                      FUN=function(x) sum(x, na.rm=T))
stats.all_FRA1$no.SNP.per.RADtag.med <- round(tapply(FRA1_pi_m$no_sites, INDEX=FRA1_pi_m$chr.type, 
                                                       FUN=function(x) median(x, na.rm=T)), 1)

stats.all_FRA1$no.SNP.per.RADtag.q1 <- round(tapply(FRA1_pi_m$no_sites, INDEX=FRA1_pi_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.25, na.rm=T)), 1)

stats.all_FRA1$no.SNP.per.RADtag.q3 <- round(tapply(FRA1_pi_m$no_sites, INDEX=FRA1_pi_m$chr.type, 
                                                      FUN=function(x) quantile(x, 0.75, na.rm=T)), 1)

stats.all_FRA1$window_length.med <- round(tapply(FRA1_pi_m$window_length, INDEX=FRA1_pi_m$chr.type, 
                                                   FUN=function(x) median(x, na.rm=T)), 1)

stats.all_FRA1$window_length.min <- round(tapply(FRA1_pi_m$window_length, INDEX=FRA1_pi_m$chr.type, 
                                                   FUN=function(x) min(x, na.rm=T)),1)

stats.all_FRA1$window_length.max <- round(tapply(FRA1_pi_m$window_length, INDEX=FRA1_pi_m$chr.type, 
                                                   FUN=function(x) max(x, na.rm=T)),1)

stats.all_FRA1$dist.med <- round(tapply(FRA1_pi_m$dist, INDEX=FRA1_pi_m$chr.type, 
                                          FUN=function(x) median(x, na.rm=T)), 0)
stats.all_FRA1$dist.q1 <- round(tapply(FRA1_pi_m$dist, INDEX=FRA1_pi_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.25, na.rm=T)),0)
stats.all_FRA1$dist.q3 <- round(tapply(FRA1_pi_m$dist, INDEX=FRA1_pi_m$chr.type, 
                                         FUN=function(x) quantile(x, 0.75, na.rm=T)),0)


write.table(stats.all_FRA1, "Stats output/stats.all.FRA1_230912.csv", sep=",")

###rolling median and rolling mean####

for(n in 1:13){
  FRA1_pi_m$movmed.51[FRA1_pi_m$chromosome==n] <- 
    rollmedian( FRA1_pi_m$avg_pi[FRA1_pi_m$chromosome==n], 51, fill=NA) 
  FRA1_pi_m$movmean.50[FRA1_pi_m$chromosome==n] <- 
    rollmean( FRA1_pi_m$avg_pi[FRA1_pi_m$chromosome==n], 50, fill=NA) 
}

str(FRA1_pi_m)



###add windowed calculations####
#non-overlapping windows
window.size <- 5000000

n <- 1
for (n in c(1:12))
{
  no <- floor(start$length[n]/window.size)
  current <- data.frame(chr = as.numeric(paste(n)), window=1:no, fst=NA, mid=NA, start=NA, snp=NA, scaff=NA, seq.len=NA)
  shift <- (start$length[n]-window.size*no)/2
  current$start <- seq(from=shift, to = shift + window.size*no-1, by=window.size)
  current$mid <- current$start + window.size/2
  chr <- FRA1_pi_m[FRA1_pi_m$chromosome==n,] 
  m <- 1
  for (m in c(1:no)){
    win <-     chr[chr$window_pos_1 >= current$start[m] 
                   & chr$window_pos_2 < current$start[m]+window.size,]
    current$fst[m] <- mean(win$avg_pi)
    current$scaff[m] <- dim(win)[1]
    current$sites[m] <- sum(win$no_sites)
    current$seq.len[m] <- sum(win$window_length)
    assign(paste("chr",n, sep="_"), current)
  }
  
}

data.wind_FRA1 <- rbind(chr_1, chr_2, chr_3, chr_4, chr_5, chr_6, chr_7, chr_8, chr_9, chr_10, chr_11, chr_12)


data.wind_FRA1$chr.type[data.wind_FRA1$chr<12] <- "Autosomes"
data.wind_FRA1$chr.type[data.wind_FRA1$chr==12] <- "X-chr."

data.wind_FRA1$chr.type <- factor(data.wind_FRA1$chr.type)


sum(data.wind_FRA1$scaff < 4, na.rm = T)
data.wind_FRA1 <- data.wind_FRA1[data.wind_FRA1$scaff > 4,]
str(data.wind_FRA1)
window.summary_FRA1 <- data.frame(chr.type=levels(data.wind_FRA1$chr.type), window.size=window.size,
                             no.RADtags.per.wind.med =NA, 
                             no.RADtags.q1 =NA, 
                             no.RADtags.q3 =NA,
                             no.sites.per.wind.med=NA,
                             no.sites.per.wind.q1=NA, 
                             no.sites.per.wind.q3=NA,
                             seq.length.per.wind.med=NA,
                             seq.length.per.wind.q1=NA,
                             seq.length.per.wind.q3=NA)

window.summary_FRA1$no.RADtags.per.wind.med <- round(tapply(data.wind_FRA1$scaff, data.wind_FRA1$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_FRA1$no.RADtags.q1 <- round(tapply(data.wind_FRA1$scaff, data.wind_FRA1$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_FRA1$no.RADtags.q3 <- round(tapply(data.wind_FRA1$scaff, data.wind_FRA1$chr.type, function(x) quantile(x, 0.75, na.rm=T)),1)

window.summary_FRA1$no.sites.per.wind.med <- round(tapply(data.wind_FRA1$sites , data.wind_FRA1$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_FRA1$no.sites.per.wind.q1 <- round(tapply(data.wind_FRA1$sites , data.wind_FRA1$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_FRA1$no.sites.per.wind.q3 <- round(tapply(data.wind_FRA1$sites , data.wind_FRA1$chr.type, function(x) quantile(x, 0.75, na.rm=T)), 1)


window.summary_FRA1$seq.length.per.wind.med <- round(tapply(data.wind_FRA1$seq.len, data.wind_FRA1$chr.type, function(x) median(x, na.rm=T)),0)
window.summary_FRA1$seq.length.per.wind.q1 <- round(tapply(data.wind_FRA1$seq.len, data.wind_FRA1$chr.type, function(x) quantile(x, 0.25,na.rm=T)),0)
window.summary_FRA1$seq.length.per.wind.q3 <- round(tapply(data.wind_FRA1$seq.len, data.wind_FRA1$chr.type, function(x) quantile(x,0.75, na.rm=T)),0)

window.summary_FRA1
write.csv(window.summary_FRA1, "Stats output/FRA1_window_stats_230908.csv")
tail(data.wind_FRA1)

###add start & end of  chromosome ####
#adjust no column as needed
str(FRA1_pi_m)
add <- matrix(nrow=26, ncol=16, NA)
add[,2] <- add.end[,3]
add[,12] <- add.end[,10]
add 
add<- data.frame(add)
colnames(add) <- colnames(FRA1_pi_m)
FRA1_pi_m <- rbind(FRA1_pi_m, add)
str(FRA1_pi_m)

table(FRA1_pi_m$chromosome)






# pi POL3####
###mid window length & mid position####
POL3_pi$window_length <- POL3_pi$window_pos_2- POL3_pi$window_pos_1
POL3_pi$mid_pos <- POL3_pi$window_pos_1 + POL3_pi$window_length/2

###remove pi < 0, pi =NA, window length < 100 and adjust chromosome and SNP var ####
table(POL3_pi$avg_pi <0)
table(POL3_pi$window_length < 100, POL3_pi$chromosome)
summary(POL3_pi)
boxplot(POL3_pi$window_length~POL3_pi$chromosome)
table(is.na(POL3_pi$avg_pi))
table(POL3_pi$avg_pi <0 | POL3_pi$window_length < 100 | is.na(POL3_pi$avg_pi)==T)

str(POL3_pi)

POL3_pi_m <- POL3_pi %>% 
  mutate(chromosome= as.numeric(gsub("chr", "", chromosome)))%>% 
  mutate(avg_pi = replace(avg_pi, which(avg_pi <0), NA)) %>% 
  mutate(avg_pi = replace(avg_pi, which(window_length < 100), NA)) %>% 
  mutate(SNP=seq(1,dim(POL3_pi)[1], by=1), .after=avg_pi) %>% 
  na.omit() 

dim(POL3_pi)[1] - dim(POL3_pi_m)[1]
str(POL3_pi_m)

#define chromosome types
POL3_pi_m$chr.type[POL3_pi_m$chromosome<12] <- "Autosomes"
POL3_pi_m$chr.type[POL3_pi_m$chromosome==12] <- "X-chr."
POL3_pi_m$chr.type[POL3_pi_m$chromosome==13] <- "Y-chr."

POL3_pi_m$chr.type <- factor(POL3_pi_m$chr.type)

###distance between SNP +stats ####

POL3_pi_m$dist <- NA

for(n in 2:dim(POL3_pi_m)[1]) 
{
  POL3_pi_m$dist[n]  <- POL3_pi_m$window_pos_1[n]-POL3_pi_m$window_pos_1[n-1]
}

sum(POL3_pi_m$dist <0, na.rm=T)
POL3_pi_m$dist[POL3_pi_m$dist<0] <- NA

summary(POL3_pi_m$dist)
stats.all_POL3 <- data.frame(comparison="POL3",chr.type=levels(POL3_pi_m$chr.type), 
                             no.RADtags.tot = NA, 
                             no.SNP.tot = NA,
                             no.SNP.per.RADtag.med=NA, 
                             no.SNP.per.RADtag.q1=NA, 
                             no.SNP.per.RADtag.q3=NA,
                             window_length.med=NA,
                             window_length.min=NA,
                             window_length.max=NA,
                             dist.med=NA,
                             dist.q1=NA,
                             dist.q3=NA
)

stats.all_POL3$no.RADtags.tot <- table(POL3_pi_m$chr.type)
stats.all_POL3$no.SNP.tot <- tapply(POL3_pi_m$no_sites, INDEX=POL3_pi_m$chr.type, 
                                    FUN=function(x) sum(x, na.rm=T))
stats.all_POL3$no.SNP.per.RADtag.med <- round(tapply(POL3_pi_m$no_sites, INDEX=POL3_pi_m$chr.type, 
                                                     FUN=function(x) median(x, na.rm=T)), 1)

stats.all_POL3$no.SNP.per.RADtag.q1 <- round(tapply(POL3_pi_m$no_sites, INDEX=POL3_pi_m$chr.type, 
                                                    FUN=function(x) quantile(x, 0.25, na.rm=T)), 1)

stats.all_POL3$no.SNP.per.RADtag.q3 <- round(tapply(POL3_pi_m$no_sites, INDEX=POL3_pi_m$chr.type, 
                                                    FUN=function(x) quantile(x, 0.75, na.rm=T)), 1)

stats.all_POL3$window_length.med <- round(tapply(POL3_pi_m$window_length, INDEX=POL3_pi_m$chr.type, 
                                                 FUN=function(x) median(x, na.rm=T)), 1)

stats.all_POL3$window_length.min <- round(tapply(POL3_pi_m$window_length, INDEX=POL3_pi_m$chr.type, 
                                                 FUN=function(x) min(x, na.rm=T)),1)

stats.all_POL3$window_length.max <- round(tapply(POL3_pi_m$window_length, INDEX=POL3_pi_m$chr.type, 
                                                 FUN=function(x) max(x, na.rm=T)),1)

stats.all_POL3$dist.med <- round(tapply(POL3_pi_m$dist, INDEX=POL3_pi_m$chr.type, 
                                        FUN=function(x) median(x, na.rm=T)), 0)
stats.all_POL3$dist.q1 <- round(tapply(POL3_pi_m$dist, INDEX=POL3_pi_m$chr.type, 
                                       FUN=function(x) quantile(x, 0.25, na.rm=T)),0)
stats.all_POL3$dist.q3 <- round(tapply(POL3_pi_m$dist, INDEX=POL3_pi_m$chr.type, 
                                       FUN=function(x) quantile(x, 0.75, na.rm=T)),0)


write.table(stats.all_POL3, "Stats output/stats.all.POL3_230912.csv", sep=",")

###rolling median and rolling mean####

for(n in 1:13){
  POL3_pi_m$movmed.51[POL3_pi_m$chromosome==n] <- 
    rollmedian( POL3_pi_m$avg_pi[POL3_pi_m$chromosome==n], 51, fill=NA) 
  POL3_pi_m$movmean.50[POL3_pi_m$chromosome==n] <- 
    rollmean( POL3_pi_m$avg_pi[POL3_pi_m$chromosome==n], 50, fill=NA) 
}

str(POL3_pi_m)



###add windowed calculations####
#non-overlapping windows
window.size <- 5000000

n <- 1
for (n in c(1:12))
{
  no <- floor(start$length[n]/window.size)
  current <- data.frame(chr = as.numeric(paste(n)), window=1:no, fst=NA, mid=NA, start=NA, snp=NA, scaff=NA, seq.len=NA)
  shift <- (start$length[n]-window.size*no)/2
  current$start <- seq(from=shift, to = shift + window.size*no-1, by=window.size)
  current$mid <- current$start + window.size/2
  chr <- POL3_pi_m[POL3_pi_m$chromosome==n,] 
  m <- 1
  for (m in c(1:no)){
    win <-     chr[chr$window_pos_1 >= current$start[m] 
                   & chr$window_pos_2 < current$start[m]+window.size,]
    current$fst[m] <- mean(win$avg_pi)
    current$scaff[m] <- dim(win)[1]
    current$sites[m] <- sum(win$no_sites)
    current$seq.len[m] <- sum(win$window_length)
    assign(paste("chr",n, sep="_"), current)
  }
  
}

data.wind_POL3 <- rbind(chr_1, chr_2, chr_3, chr_4, chr_5, chr_6, chr_7, chr_8, chr_9, chr_10, chr_11, chr_12)


data.wind_POL3$chr.type[data.wind_POL3$chr<12] <- "Autosomes"
data.wind_POL3$chr.type[data.wind_POL3$chr==12] <- "X-chr."

data.wind_POL3$chr.type <- factor(data.wind_POL3$chr.type)


sum(data.wind_POL3$scaff < 4, na.rm = T)
data.wind_POL3 <- data.wind_POL3[data.wind_POL3$scaff > 4,]
str(data.wind_POL3)
window.summary_POL3 <- data.frame(chr.type=levels(data.wind_POL3$chr.type), window.size=window.size,
                                  no.RADtags.per.wind.med =NA, 
                                  no.RADtags.q1 =NA, 
                                  no.RADtags.q3 =NA,
                                  no.sites.per.wind.med=NA,
                                  no.sites.per.wind.q1=NA, 
                                  no.sites.per.wind.q3=NA,
                                  seq.length.per.wind.med=NA,
                                  seq.length.per.wind.q1=NA,
                                  seq.length.per.wind.q3=NA)

window.summary_POL3$no.RADtags.per.wind.med <- round(tapply(data.wind_POL3$scaff, data.wind_POL3$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_POL3$no.RADtags.q1 <- round(tapply(data.wind_POL3$scaff, data.wind_POL3$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_POL3$no.RADtags.q3 <- round(tapply(data.wind_POL3$scaff, data.wind_POL3$chr.type, function(x) quantile(x, 0.75, na.rm=T)),1)

window.summary_POL3$no.sites.per.wind.med <- round(tapply(data.wind_POL3$sites , data.wind_POL3$chr.type, function(x) median(x, na.rm=T)), 1)
window.summary_POL3$no.sites.per.wind.q1 <- round(tapply(data.wind_POL3$sites , data.wind_POL3$chr.type, function(x) quantile(x, 0.25, na.rm=T)), 1)
window.summary_POL3$no.sites.per.wind.q3 <- round(tapply(data.wind_POL3$sites , data.wind_POL3$chr.type, function(x) quantile(x, 0.75, na.rm=T)), 1)


window.summary_POL3$seq.length.per.wind.med <- round(tapply(data.wind_POL3$seq.len, data.wind_POL3$chr.type, function(x) median(x, na.rm=T)),0)
window.summary_POL3$seq.length.per.wind.q1 <- round(tapply(data.wind_POL3$seq.len, data.wind_POL3$chr.type, function(x) quantile(x, 0.25,na.rm=T)),0)
window.summary_POL3$seq.length.per.wind.q3 <- round(tapply(data.wind_POL3$seq.len, data.wind_POL3$chr.type, function(x) quantile(x,0.75, na.rm=T)),0)

window.summary_POL3
write.csv(window.summary_POL3, "Stats output/POL3_window_stats_230912.csv")
tail(data.wind_POL3)

###add start & end of  chromosome ####
#adjust no column as needed
str(POL3_pi_m)
add <- matrix(nrow=26, ncol=16, NA)
add[,2] <- add.end[,3]
add[,12] <- add.end[,10]
add 
add<- data.frame(add)
colnames(add) <- colnames(POL3_pi_m)
POL3_pi_m <- rbind(POL3_pi_m, add)
str(POL3_pi_m)

table(POL3_pi_m$chromosome)










# FIGURE Fig. S14####--------------------------------
#example for individual plot:
  #pdf("FRA1_pi_window.pdf", height=4, width=10)
  #jpeg("FRA1_pi_window.jpg", height=4, width=10, units="in", res=300)
  #dev.off()

jpeg("Figures/Fig_S14_241210.jpg", height=8, width=6.5, units="in", res=400)

#multi-panel plot pars
par(mfrow=c(5,1), xpd=T, mar=c(4, 5, 2, 2), oma=c(1,0,1,0))

#### (a) pi FRA1 ####
str(FRA1_pi_m)
table(FRA1_pi_m$chromosome)
manhattan(FRA1_pi_m, 
          chr="chromosome",bp="mid_pos",p = "avg_pi", 
          snp = "SNP", logp=FALSE, xlab="", 
          ylab=expression(paste(pi)), 
          col = alpha(c("blue","black"), 0.15), 
          annotateTop = F, 
          suggestiveline= F,
          cex=0.1, xaxt="n", yaxt="n",
          bty="n", ylim=c(-0.001,0.05), cex.lab=1.5, xpd=F)


axis(side=1, labels=c(1:11, "X", "Y"), at=start$tick, tick=F, line=F, pos=0)
axis(side=2, labels=seq(from=0.0, to=0.05, by=0.01), at=seq(from=0.0, to=0.05, by=0.01), las=1)
mtext(expression(italic("S. latifolia")~(FRA1)), cex = 0.9, line=0.5, adj=0, font.main=1, at=0)
addfiglab("A")
str(start)

for (n in c(1,3,5,7,9,11)){
  lines( sum(start$length[start$chromosome <n]) + data.wind_FRA1$mid[data.wind_FRA1$chr==n], data.wind_FRA1$fst[data.wind_FRA1$chr==n]
         , col=c("blue"), lwd=1.5)
}

str(data.wind_FRA1)
for (n in c(2,4,6,8,10,12)){
  lines( sum(start$length[start$chromosome <n]) + data.wind_FRA1$mid[data.wind_FRA1$chr==n], data.wind_FRA1$fst[data.wind_FRA1$chr==n]
         , col=c("black"), lwd=1.5)
}

points(sum(start$length[start$chromosome <13])+FRA1_pi_m$mid_pos[FRA1_pi_m$chromosome==13], 
       FRA1_pi_m$avg_pi[FRA1_pi_m$chromosome==13], 
       col="blue", pch=21, bg="blue", cex=0.4)


for (n in c(1,3,5,7,9,11,13)){
  points(sum(start$length[start$chromosome <n])+FRA1_pi_m$mid_pos[FRA1_pi_m$avg_pi>= 0.05 & FRA1_pi_m$chromosome==n], 
         rep(0.05,time=length(FRA1_pi_m$mid_pos[FRA1_pi_m$avg_pi>= 0.05 & FRA1_pi_m$chromosome==n])), 
         col=alpha(c("blue"), 0.3), pch=21, bg=alpha(c("blue"), 0.3), cex=0.1)
}


for (n in c(2,4,6,8,10,12)){
  points(sum(start$length[start$chromosome <n])+FRA1_pi_m$mid_pos[FRA1_pi_m$avg_pi>= 0.05 & FRA1_pi_m$chromosome==n], 
         rep(0.05,time=length(FRA1_pi_m$mid_pos[FRA1_pi_m$avg_pi>= 0.05 & FRA1_pi_m$chromosome==n])), 
         col=alpha(c("black"), 0.3), pch=21, bg=alpha(c("black"), 0.3), cex=0.1)
}

#### (b) pi POL3  ####
str(POL3_pi_m)
summary(POL3_pi_m)
sum(POL3_pi_m$avg_pi>0.05, na.rm=T)
table(POL3_pi_m$avg_pi>0.05, POL3_pi_m$chromosome)
manhattan(POL3_pi_m, 
          chr="chromosome",bp="mid_pos",p = "avg_pi", 
          snp = "SNP", logp=FALSE, xlab="", 
          ylab=expression(paste(pi)), 
          col = alpha(c("blue","black"), 0.15), 
          annotateTop = F, 
          suggestiveline= F,
          cex=0.1, xaxt="n", yaxt="n", xpd=F,
          bty="n", ylim=c(-0.001,0.05), cex.lab=1.5)


axis(side=1, labels=c(1:11, "X", "Y"), at=start$tick, tick=F, line=F, pos=0)
axis(side=2, labels=seq(from=0.0, to=0.05, by=0.01), at=seq(from=0.0, to=0.05, by=0.01), las=1)

mtext(expression(italic("S. latifolia")~(POL3)), cex = 0.9, line= 0.5, adj=0, font.main=1, at=0)
addfiglab("B")
str(data.wind_POL3)

for (n in c(1,3,5,7,9,11)){
  lines( sum(start$length[start$chromosome <n]) + data.wind_POL3$mid[data.wind_POL3$chr==n], data.wind_POL3$fst[data.wind_POL3$chr==n]
         , col=c("blue"), lwd=1.5)
}
for (n in c(2,4,6,8,10,12)){
  lines( sum(start$length[start$chromosome <n]) + data.wind_POL3$mid[data.wind_POL3$chr==n], data.wind_POL3$fst[data.wind_POL3$chr==n]
         , col=c("black"), lwd=1.5)
}

points(sum(start$length[start$chromosome <13])+POL3_pi_m$mid_pos[POL3_pi_m$chromosome==13], 
       POL3_pi_m$avg_pi[POL3_pi_m$chromosome==13], 
       col="blue", pch=21, bg="blue", cex=0.4)


for (n in c(1,3,5,7,9,11,13)){
  points(sum(start$length[start$chromosome <n])+POL3_pi_m$mid_pos[POL3_pi_m$avg_pi>= 0.05 & POL3_pi_m$chromosome==n], 
         rep(0.05,time=length(POL3_pi_m$mid_pos[POL3_pi_m$avg_pi>= 0.05 & POL3_pi_m$chromosome==n])), 
         col=alpha(c("blue"), 0.3), pch=21, bg=alpha(c("blue"), 0.3), cex=0.1)
}


for (n in c(2,4,6,8,10,12)){
  points(sum(start$length[start$chromosome <n])+POL3_pi_m$mid_pos[POL3_pi_m$avg_pi>= 0.05 & POL3_pi_m$chromosome==n], 
         rep(0.05,time=length(POL3_pi_m$mid_pos[POL3_pi_m$avg_pi>= 0.05 & POL3_pi_m$chromosome==n])), 
         col=alpha(c("black"), 0.3), pch=21, bg=alpha(c("black"), 0.3), cex=0.1)
}

#### (c) Fst within S. latifolia (FRAPOL)####
str(FRA1POL3_fst_m.x)
table(FRA1POL3_fst_m.x$chromosome)
FRA1POL3_fst_m.x$chromosome <- as.numeric(FRA1POL3_fst_m.x$chromosome)
manhattan(FRA1POL3_fst_m.x, 
          chr="chromosome",bp="mid_pos",p = "avg_hudson_fst", 
          snp = "SNP", logp=FALSE, xlab="", 
          ylab=expression("F"[ST]), 
          col = alpha(c("blue","black"), 0.15), 
          annotateTop = F, 
          suggestiveline= F,
          cex=0.2, xaxt="n", 
          bty="n", ylim=c(-0.05,1.1), cex.lab=1.3)

axis(side=1, labels=c(1:11, "X", "Y"), at=start$tick, tick=F, line=F, pos=0)
mtext(expression(within~italic("S. latifolia")~(FRA1-POL3)), cex = 0.9, line=0, adj=0, font.main=1, at=0)
addfiglab("C")


for (n in c(1,3,5,7,9,11,13)){
  lines( sum(start$length[start$chr <n])+data.wind_FRAPOL$mid[data.wind_FRAPOL$chr==n], data.wind_FRAPOL$fst[data.wind_FRAPOL$chr==n]
         , col=c("blue"), lwd=1.5)
}


for (n in c(2,4,6,8,10,12)){
  lines( sum(start$length[start$chr <n]) + data.wind_FRAPOL$mid[data.wind_FRAPOL$chr==n], data.wind_FRAPOL$fst[data.wind_FRAPOL$chr==n]
         , col=c("black"), lwd=1.5)
}

points(sum(start$length[start$chr <13])+FRA1POL3_fst_m.x$mid_pos[FRA1POL3_fst_m.x$chromosome==13], 
       FRA1POL3_fst_m.x$avg_hudson_fst[FRA1POL3_fst_m.x$chromosome==13], 
       col="blue", pch=21, bg="blue", cex=0.4)

#### (d) Fst btw BELFRA ####
str(BEL1FRA1_fst_m.x)
table(BEL1FRA1_fst_m.x$chromosome)
BEL1FRA1_fst_m.x$chromosome <- as.numeric(BEL1FRA1_fst_m.x$chromosome)
manhattan(BEL1FRA1_fst_m.x, 
          chr="chromosome",bp="mid_pos",p = "avg_hudson_fst", 
          snp = "SNP", logp=FALSE, xlab="", 
          ylab=expression("F"[ST]), 
          col = alpha(c("blue","black"), 0.2), 
          annotateTop = F, 
          suggestiveline= F,
          cex=0.2, xaxt="n", #main=" 5 Mb between species FST (SL-SD)", 
          bty="n", ylim=c(-0.05,1.1), cex.lab=1.3)

axis(side=1, labels=c(1:11, "X", "Y"), at=start$tick, tick=F, line=F, pos=0)


mtext(expression(between~italic("S. latifolia")~(FRA1)~and~italic("S. dioica")~(BEL1)), cex = 0.9, line=0, adj=0, font.main=1, at=0)
addfiglab("D")

for (n in c(1,3,5,7,9,11)){
  lines( sum(start$length[start$chr <n]) + data.wind_BELFRA$mid[data.wind_BELFRA$chr==n], data.wind_BELFRA$fst[data.wind_BELFRA$chr==n]
         , col=c("blue"), lwd=1.5)
}
for (n in c(2,4,6,8,10,12)){
  lines( sum(start$length[start$chr <n]) + data.wind_BELFRA$mid[data.wind_BELFRA$chr==n], data.wind_BELFRA$fst[data.wind_BELFRA$chr==n]
         , col=c("black"), lwd=1.5)
}

points(sum(start$length[start$chr <13])+BEL1FRA1_fst_m.x$mid_pos[BEL1FRA1_fst_m.x$chromosome==13], 
       BEL1FRA1_fst_m.x$avg_hudson_fst[BEL1FRA1_fst_m.x$chromosome==13], 
       col="blue", pch=21, bg="blue", cex=0.4)



#### (e) Fst btw BELPOL ####
BEL1POL3_fst_m.x$chromosome <- as.numeric(BEL1POL3_fst_m.x$chromosome)
manhattan(BEL1POL3_fst_m.x, 
          chr="chromosome",bp="mid_pos",p = "avg_hudson_fst", 
          snp = "SNP", logp=FALSE, xlab="", 
          ylab=expression("F"[ST]), 
          col = alpha(c("blue","black"), 0.2), 
          annotateTop = F, 
          suggestiveline= F,
          cex=0.2, xaxt="n", #main="between species FST (SL-SD)", 
          bty="n", ylim=c(-0.05,1.1), cex.lab=1.3)

axis(side=1, labels=c(1:11, "X", "Y"), at=start$tick, tick=F, line=F, pos=0, cex.lab=1.3)
mtext(expression(between~italic("S. latifolia")~(POL3)~and~italic("S. dioica")~(BEL1)), cex = 0.9, line=0, adj=0, font.main=1, at=0)
addfiglab("E")

for (n in c(1,3,5,7,9,11)){
  lines( sum(start$length[start$chr <n]) + data.wind_BELPOL$mid[data.wind_BELPOL$chr==n], data.wind_BELPOL$fst[data.wind_BELPOL$chr==n]
         , col=c("blue"), lwd=1.5)
}
for (n in c(2,4,6,8,10,12)){
  lines( sum(start$length[start$chr <n]) + data.wind_BELPOL$mid[data.wind_BELPOL$chr==n], data.wind_BELPOL$fst[data.wind_BELPOL$chr==n]
         , col=c("black"), lwd=1.5)
}

points(sum(start$length[start$chromosome <13])+BEL1POL3_fst_m.x$mid_pos[BEL1POL3_fst_m.x$chromosome==13], 
       BEL1POL3_fst_m.x$avg_hudson_fst[BEL1POL3_fst_m.x$chromosome==13], 
       col="blue", pch=21, bg="blue", cex=0.4)

####overall x axis label####
mtext("Chromosome", side=1, line=-1, cex=1, outer=T, adj=NA, font=2)

dev.off()



# window data for export####
write.table(data.wind_BELFRA, "data/wind_5Mb_BELFRA_Fst_230908.csv", sep=",")
write.table(data.wind_BELPOL, "data/wind_5Mb_BELPOL_Fst.csv", sep=",")
write.table(data.wind_FRAPOL, "data/wind_5Mb_FRAPOL_Fst.csv", sep=",")
write.table(data.wind_FRA1, "data/wind_5Mb_FRA1_pi_230908.csv", sep=",")
write.table(data.wind_POL3, "data/wind_5Mb_POL3_pi.csv", sep=",")


str(BEL1FRA1_fst_m)

write.table(BEL1FRA1_fst_m[BEL1FRA1_fst_m$chromosome=="13", 1:11],"data/BEL1FRA1_fst_Y_230908.csv", sep="," )
str(FRA1_pi_m)

write.table(FRA1_pi_m[FRA1_pi_m$chromosome=="13", 1:13],"data/FRA1_pi_Y_230908.csv", sep="," )

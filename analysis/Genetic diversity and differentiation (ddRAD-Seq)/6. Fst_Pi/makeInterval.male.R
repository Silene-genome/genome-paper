#!/usr/bin/env Rscript
# Rscript ../script/makeInterval.R /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/female/chr2_2pop_noINDEL_minDP6_geno05perPop.DP.tsv 2

args = commandArgs(trailingOnly=TRUE)

infile= args[1]
chr = args[2]
outfile=args[3]
t=read.table(infile,sep="\t",header=F)
table = t
pos = as.numeric(table[,2])

interval = split(pos, cumsum(c(1, diff(pos))!=1))

intervals = matrix(NA,ncol=3,nrow=length(interval))
ii = 1
for (i in names(interval)) {
    intervals[ii,1] = ii
    pos_v = interval[[i]]
    if (length(pos_v)<=1) {
            intervals[ii,2] =pos_v
            intervals[ii,3] = pos_v
    } else {
            intervals[ii,2] = pos_v[1]
            intervals[ii,3] = pos_v[length(pos_v)]
    }
    ii = ii+1
}

t<-intervals
span = 125

windowstart = c()
windowend = c()
bases = c()
for(i in 1:nrow(t)) {
  row.start = t[i,2]
  row.end = t[i,3]
  if (i==1) {
          start=t[i,2]
          Nbases = row.end - row.start + 1
  }

  if(( row.start > (start+span) )&&(row.end > (start+span)) ) { # out of old range                                                                                                                                                                                                                                                                                      
    #cat(start,"\t",start+span,"\t",Nbases,"\n")                                                                                                                                                                                                                                                                                                                        
    windowstart = append(windowstart,start)
    windowend = append(windowend,t[i-1,3])
    bases = append(bases,Nbases)
    start = row.start
    Nbases = row.end - row.start + 1

  } else if ( (row.start < (start+span) ) && (row.end > (start+span)) ) { # part of end of range fall into row interval                                                                                                                                                                                                                                                 
    cat("Warning: ", row.start, " range ends in the row interval","\n")
    rowbases = row.end - row.start+1
    Nbases = Nbases + rowbases
  } else if ( row.start == start) {
    next
  } else { # interval (row) falls into the range                                                                                                                                                                                                                                                                                                                        
    rowbases = row.end - row.start+1
    Nbases = Nbases + rowbases
  }

  if (i==nrow(t)) { #end of t                                                                                                                                                                                                                                                                                                                                           
     #cat(start,"\t",start+span,"\t",Nbases,"\n")                                                                                                                                                                                                                                                                                                                       
     windowstart = append(windowstart,start)
     windowend =append(windowend,row.end)
     bases = append(bases,Nbases)


  }
}

dynamicWindow = data.frame(chr=chr,start=windowstart,end=windowend+1,N=bases)

thres_interval = 100
write.table(dynamicWindow[which(dynamicWindow$N>=thres_interval),],file=outfile,sep="\t",col.names=F,row.names=F,quote=FALSE)
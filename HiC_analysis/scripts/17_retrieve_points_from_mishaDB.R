#!/usr/bin/Rscript

require(data.table)
require(misha)
require(plyr)
require(dplyr)
require(reshape2)

options(max.print=500000000, gmax.data.size=100000000, gmax.mem.usage=100000000)
#options()

### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()
print(dbDir)

options(scipen=12, gmax.data.size=5e7)

track <- 'XXXtrackXXX'

### Load chromosome sizes
chromSizes	<- read.table(paste0(dbDir,'/chrom_sizes.txt'), header=F)
colnames(chromSizes) <- c('chr','len')
#print(chromSizes)

### Define the plot parameters
nchunk <- XXXnchunkXXX

chrom1      <- "XXXchrom1XXX"
chromSize  <- dplyr::filter(chromSizes,chromSizes$chr == gsub('chr','',chrom1))$len #; print(chromSize)
start1      <- XXXstart1XXX
if(start1 < 0){start1 = 0}
end1        <- XXXend1XXX
if(end1 > chromSize){end1 = chromSize}
print(paste0(chrom1," ",start1," ",end1))  

chrom2      <- "XXXchrom2XXX"
chromSize  <- dplyr::filter(chromSizes,chromSizes$chr == gsub('chr','',chrom2))$len #; print(chromSize)
start2      <- XXXstart2XXX
if(start2 < 0){start2 = 0}
end2        <- XXXend2XXX
if(end2 > chromSize){end2 = chromSize}
print(paste0(chrom2," ",start2," ",end2))

outfile <- paste0("Observed_merged_",chrom1,"_",start1,"_",end1,"bp_",chrom2,"_",start2,"_",end2,"bp_chunk_",nchunk,".ints")
print(outfile)

### Define a 2d interval
interval <- gintervals.2d(chrom1,start1,end1,chrom2,start2,end2)

scores   <- gextract(track,interval,colnames='score') #,band=c(-(end-start),-10e3))
scores[is.na(scores)] <- 0

### Write output
write.table(scores, outfile, append = FALSE, sep = " ", dec = ".",row.names = F, col.names = T, quote=FALSE)
quit()

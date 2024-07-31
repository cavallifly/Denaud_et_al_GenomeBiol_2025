library(misha)
library(shaman)

source("./scripts/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

### misha working DB
mDBloc <-  "./mishaDB/trackdb/"
db     <- "dm6"
dbDir  <- paste0(mDBloc,db,"/")
gdb.init(dbDir)
gdb.reload()

tssCoordinates           <- gintervals.load("intervals.ucscCanTSS")
rownames(tssCoordinates) <- tssCoordinates$geneName

insTracks <- list(
        larvae_DWT_100    = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w100kb_r2kb',
        larvae_DWT_150    = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w150kb_r2kb',
        larvae_DWT_200    = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w200kb_r2kb',
        larvae_DWT_250    = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w250kb_r2kb',
        larvae_DWT_300    = 'insulation.INS_hic_larvae_DWT_merge_dm6_BeS_w300kb_r2kb',	
	larvae_DPRE1_100 = 'insulation.INS_hic_larvae_DPRE1_merge_dm6_BeS_w100kb_r2kb',
	larvae_DPRE1_150 = 'insulation.INS_hic_larvae_DPRE1_merge_dm6_BeS_w150kb_r2kb',
	larvae_DPRE1_200 = 'insulation.INS_hic_larvae_DPRE1_merge_dm6_BeS_w200kb_r2kb',
	larvae_DPRE1_250 = 'insulation.INS_hic_larvae_DPRE1_merge_dm6_BeS_w250kb_r2kb',
	larvae_DPRE1_300 = 'insulation.INS_hic_larvae_DPRE1_merge_dm6_BeS_w300kb_r2kb',
	larvae_DPRE1Up_100 = 'insulation.INS_hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS_w100kb_r2kb',
	larvae_DPRE1Up_150 = 'insulation.INS_hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS_w150kb_r2kb',
	larvae_DPRE1Up_200 = 'insulation.INS_hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS_w200kb_r2kb',
	larvae_DPRE1Up_250 = 'insulation.INS_hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS_w250kb_r2kb',
	larvae_DPRE1Up_300 = 'insulation.INS_hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS_w300kb_r2kb'	
        )
refSet = "larvae_DWT"

# Genomic locations of interest
pre1     = gintervals("chr2L",16422514,16422515)
boundary1 <- gintervals("chr2L",16354000,16354001)
boundary2 <- gintervals("chr2L",16494000,16494001)
dacTSS   = tssCoordinates["dac",1:3]

chr <- "chr2L"
center <- 16485998
insIntrv <- gintervals(chr, center - 165000, center + 25000)

# Retrieve the insulation tracks from mishaDB
print("Get the min, mean, stddev, and max of the insulation tracks")
insData <- list()
for(set in names(insTracks))
{
    print(set, row.names=F, quote=F)
    track <- insTracks[[set]]
    data <- gextract(track,insIntrv,iterator=3e3,colnames="INS")    
    data$INS  <- clipQuantile(-data$INS,0.999)
    data$INS  <- scaleData(data$INS,0,1,1e-9)
    
    trackName  <- gsub("_hic","",gsub("_merge_dm6DPRE1Up_BeS","",gsub("_merge_dm6_BeS","",gsub("insulation.","",set))))
    trackName1 <- unlist(strsplit(trackName,"_"))[1:2]
    trackName1 <- paste(trackName1,collapse="_")
    window     <- paste0("w",unlist(strsplit(trackName,"_"))[3])
    if(length(insData[[trackName1]]) == 0)
    {
        insData[[trackName1]] <- data[c(1,2,3,4)]
	colnames(insData[[trackName1]]) <- c("chrom","start","end",window)
    } else {
        cNames <- c(colnames(insData[[trackName1]]),window)
        insData[[trackName1]] <- cbind(insData[[trackName1]],data$INS)
	colnames(insData[[trackName1]]) <- cNames
    }

}

print("")
for(trackName in names(insData))
{
    print(trackName)
    print(head(insData[[trackName]]))

    insData[[trackName]]$min    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,min)    
    insData[[trackName]]$mean   <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,mean)
    insData[[trackName]]$stddev <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,sd)
    insData[[trackName]]$max    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,max)

    outProfile <- paste0("insulationProfiles_Fig2D_",trackName,".tsv")
    write.table(insData[[trackName]], file = outProfile, sep="\t", row.names=FALSE, quote=FALSE)

    print("", row.names=F, quote=F)
}
#print(head(insData))

# Obtain the location on the plot of the GLoI
dacPos       <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - dacTSS$start))/(insIntrv$end - insIntrv$start))
pre1Pos      <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - pre1$start))/(insIntrv$end - insIntrv$start))
boundary1Pos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - boundary1$start))/(insIntrv$end - insIntrv$start))
boundary2Pos <- nrow(insData[[1]]) * (((insIntrv$end - insIntrv$start) - (insIntrv$end - boundary2$start))/(insIntrv$end - insIntrv$start))

cols     <- rainbow(length(names(insData)))
fillCols <- rainbow(length(names(insData)))
# larvae_DWT
cols[[1]]     <- rgb(  0,   0,   0, max = 255, alpha = 255)
fillCols[[1]] <- rgb(  0,   0,   0, max = 255, alpha =  50)
# DPRE1
cols[[2]]     <- rgb(  0, 255,   0, max = 255, alpha = 255)
fillCols[[2]] <- rgb(  0, 255,   0, max = 255, alpha =  50)
# DPRE1Up
cols[[3]]     <- rgb(  0, 255, 255, max = 255, alpha = 255)
fillCols[[3]] <- rgb(  0, 255, 255, max = 255, alpha =  50)
print(cols) ; print(fillCols)

print(cols) ; print(fillCols)

for(set in names(insData))
{

    if(set == "larvae_DWT")
    {
        next
    }

    if(set == "larvae_DPRE1")
    {
        pdf(paste0("insulationProfiles_Fig2D_",set,".pdf"),width=19,height=7)    
	plot(insData[[1]]$mean, type="l", ylim=c(0,1),xlim=c(0,nrow(insData[[1]])+nrow(insData[[1]])*0.15),xlab="", ylab="INSULATION", xaxt="n", frame=F)
	polygon(c(1:nrow(insData[[refSet]]), nrow(insData[[refSet]]):1), c(insData[[refSet]]$mean-insData[[refSet]]$stddev,rev(insData[[refSet]]$mean+insData[[refSet]]$stddev)), col = fillCols[which(names(insData) == refSet)], border=F)

	axis(1,c(insIntrv$start,"PRE1","TAD1",insIntrv$end),at=c(1,pre1Pos,boundary1Pos,nrow(insData[[1]])))	
        legend("topright",legend=names(insData[c(1,2)]), col=cols[c(1,2)], pch=19)
    }

    if(set == "larvae_DPRE1Up")
    {
        pdf(paste0("insulationProfiles_Fig2D_",set,".pdf"),width=19,height=7)    
	plot(insData[[1]]$mean, type="l", ylim=c(0,1),xlim=c(0,nrow(insData[[1]])+nrow(insData[[1]])*0.15),xlab="", ylab="INSULATION", xaxt="n", frame=F)
	polygon(c(1:nrow(insData[[refSet]]), nrow(insData[[refSet]]):1), c(insData[[refSet]]$mean-insData[[refSet]]$stddev,rev(insData[[refSet]]$mean+insData[[refSet]]$stddev)), col = fillCols[which(names(insData) == refSet)], border=F)

	axis(1,c(insIntrv$start,"PRE1","TAD1",insIntrv$end),at=c(1,pre1Pos,boundary1Pos,nrow(insData[[1]])))	
        legend("topright",legend=names(insData[c(1,2)]), col=cols[c(1,2)], pch=19)
    }


    print(set)

    layout(matrix(1:1,ncol=1,byrow=TRUE),respect=FALSE)
    par(mai=rep(1.5,4), cex=1.5)

    abline(v=dacPos, col=4, lwd=3, lty=2)
    abline(v=pre1Pos, col=8, lwd=3, lty=2)
    abline(v=boundary1Pos, col=5, lwd=3, lty=2)
    abline(v=boundary2Pos, col=5, lwd=3, lty=2)

    lines(insData[[set]]$mean, type="l", col=cols[which(names(insData) == set)], lwd=3)
    polygon(c(1:nrow(insData[[set]]), nrow(insData[[set]]):1), c(insData[[set]]$mean-insData[[set]]$stddev,rev(insData[[set]]$mean+insData[[set]]$stddev)), col = fillCols[which(names(insData) == set)], border=F)

    dev.off()    
}

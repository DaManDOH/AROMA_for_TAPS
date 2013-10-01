## Tumor Aberration Prediction Suite.

## This file contains all functions used to run (and used by) Tumor Aberration Prediction Suite (TAPS).
## Author: Markus Rasmussen @ Uppsala University. (Developed September 2009 --> January 2011)
## Contact: markus.rasmussen@medsci.uu.se or arrayplatform@medsci.uu.se
## Updates: Check http://array.medsci.uu.se/TAPS as we intend to extend some functionality.


## -----Instructions, allele-specific copy number analysis:--------------------------------

## Note to users:
## TAPS is a tool to help you actively investigate genomic aberrations of the most complex tumor samples.
## TAPS visualizes samples using log-ratio and allelic imbalance ratio. If the sample is deemed suitable, it can be subjected to automated copy number calling.
## Data: Validated automatic copy number calling is available for Affymetrix SNP6 and 250k/500k data.  Automatic copy number calling is suitable but unvalidated for OncoScan (an Affymetrix SNP array)

## Data preparation.
## Nexus (BioDiscovery) 5.0: find your Samples folder after segmentation (SNP-rank) and normalization.
## Other: See example package folders. Store Log-ratio as probes.txt and allele frequency as snps.txt
## in one folder per sample. Do not store other folders than sample folders in your "Samples" folder.
## The 'Array' column in these files is uninformative and may be omitted.
## Optionally (recommended) also put segments.txt in your sample folder. Please inspect the examples
## for table structure and headers.

## Prerequisites.
## Install ggplot2 from CRAN.
## If your data is not segmented already, install DNACopy (Bioconductor).

## Workflow.
## 1. Source this file.
## 2. From the folder containing your samples (sample folders) run AROMA_TAPS_plot().
## 3. Investigate the scatter plots generated in your sample folders
## 4. Investigate "<sample>__TAPS_summary.png". Do segment clusters appear to be labeled with the correct copy number?
## 5. If incorrect, an improved interpretation of Log-Ratio @ copy number 2, and the difference in
##    Log-Ratio to a deletion or amplification (+1/-1 copy) may be supplied in sampleInfo_TAPS.txt.
## 6. Run TAPS_call().
## 7. Check "<sample>__TAPS_summary.png" again if you wish, and check the new chromosome-wise images.
## 8. If all looks reasonable, you will find good copy number estimates in '<sample>_segmentCN_merged.txt'
## 9. Be wary of the result on sex chromosomes which may be difficult to auto-interpret.
## 10. Watch all images for signs of segmentation failure and tumor cell heterogeneity.


## How to use wrapped AROMA CRMAV2 to work quickly from .CEL files:
## AROMA must be installed, setup must be completed correctly, please see
## http://www.aroma-project.org/getstarted.
## Note that CEL files should reside in the AROMA subfolder rawData/[projectName]/[arrayType].
## Run AROMA_CRMAV2_toTAPS() from within the folder with the CEL files (or specify the folder as argument 1).
## You may also specify FALSE as argument 2 to NOT use the average file as a reference (if very few samples)
## Subfolders will be made with TAPS-compatible data. Once this function is finished,
## go ahead with AROMA_TAPS_plot().

library('stats')
library('fields')
AROMA_TAPS_plot <- function(directory=getwd(),xlim=c(-1,2),ylim=c(0,1),minseg=1,maxCn=8,bin=400) {
  ## This function takes a directory as input, then builds short-segment TAPS scatter plots for each sample (subdirectory) in the directory.
  ## In addition, the Log-R of copy number 2 and the difference in Log-R to deletion are estimated from the data and saved in sampleInfo_TAPS.txt.
  ## Before proceeding with TAPS_call, investigate "<sample>__TAPS_summary.png", on which TAPS interpretation of the scatter plots is marked. 
  ## If TAPS appears to be in error, you need to look carefully at the scatter plots and estimate the Log-R of copy number 2 and the difference 
  ## in Log-R to (autosome) deletion, then update sampleInfo_TAPS.txt.
  setwd(directory)
  subs <- getSubdirs()
  if (is.null(subs)) {			 							## check samples = subdirectories or a single sample = current directory
    subs=thisSubdir()
    setwd('..')
  }
  for (i in 1:length(subs)) {
    setwd(subs[i])
    name <- subs[i]
    cat(' ..loading', subs[i])
    Log2 <- AROMA_readLog2() 										## Log-R
    
    alf <- readAlf() 										## Allele Frequency
    segments <- readSegments() 								## segments if available (CBS recommended)
    
    segments$Value <- segments$Value-median(Log2$Value) 	## Median-centering
    Log2$Value <- Log2$Value-median(Log2$Value) 			## Median-centering
    
    cat(' ..processing')
    if (is.null(segments)) { 								## segmentation using DNA-copy if needed (must then be installed)
      segments <- AROMA_segment_DNAcopy(Log2)
      save.txt(segments,'_segments.txt') 
    }
    
    allRegions <- makeRegions(Log2, alf, segments)			## Calculates necessary data for segments (all functions are in this file)
    save(allRegions,file='allRegions.Rdata')
    regs <- AROMA_regsFromSegs(Log2,alf,segments,bin=bin,min=5)	## Calculates the same data for shortened segments, very useful for visualization
    
    #AROMA_scatterPlot_sample(allRegions$regions,xlim=xlim,ylim=ylim,name=name,minseg=minseg,cn=F) ## a scatter plot of the whole sample
    
    sampleInfo <- NULL
    #try( sampleInfo <- load.txt('sampleInfo_TAPS.txt'), silent=T ) ## if sampleInfo is already in the folder, no estimate is made.
    cat('..estimating ')
    if (is.null(sampleInfo)) {
      try( sampleInfo <- AROMA_getEstimates(Log2,alf,allRegions$regions), silent=T) ## Estimates the Log-R of copy number 2 and the difference to a deletion or amplification
    }
    if (is.null(sampleInfo)) { 
      sampleInfo <- data.frame('cn2'=0,'delta'=0.4)			## Default if the estimation fails
      #save.txt(sampleInfo,'sampleInfo_TAPS.txt')
      #break
    }
    save.txt(sampleInfo,'sampleInfo_TAPS.txt')
    
    cat('copy numbers')
    t <- NULL
    try (t <- AROMA_findCNs(Log2,alf,allRegions,dmin=0.9,maxCn=maxCn,ceiling=1,shift=sampleInfo$cn2,delta=sampleInfo$delta), silent=T) ## estimates the Log-R and Allelic Imbalance Ration of all variants up to maxCn
    
    try ( if (!is.null(t)) {
      u <- AROMA_setCNs(allRegions,t$int,t$ai,maxCn)				## Assigns copy number variant for all segments
      allRegions$regions <- u$regions				
      cn <- T
    }else cn <- F , silent=T)
    
    cat('..plotting.\n')
    AROMA_scatterPlot_sample(allRegions$regions,t$int,t$ai,xlim=xlim,ylim=ylim,name=name,minseg=minseg,cn=cn) ## whole-sample plot visualizing TAPS' estimates.
    for (chr in 1:23) AROMA_chromPlot_short(chr, Log2, alf, regs, xlim,ylim,name)				## Chromosome-wise plots for manual analysis
    cat('..done\n')
    setwd('..')
  }
}
###
AROMA_TAPS_call <- function(directory=getwd(),xlim=c(-1,1),ylim=c(0,1),minseg=0.2,maxCn=10) {
  ## AROMA_TAPS_call outputs the total and minor allele copy numbers of all segments as a text file, and as images for visual confirmation.
  ## sampleInfo_TAPS.txt must be present in each sample folder. If TAPS_plot could not make a good guess of the Log-R of copy number 2 
  ## and the Log-R difference to a deletion, you must interpret the scatter plots and edit sampleInfo_TAPS.txt.
  setwd(directory)
  subs <- getSubdirs()
  if (is.null(subs)) {
    subs=thisSubdir()
    setwd('..')
  }
  for (i in 1:length(subs)) {
    setwd(subs[i])
    name <- subs[i]
    cat(' ..loading', subs[i])
    Log2 <- AROMA_readLog2()
    alf <- readAlf(localDir)
    segments <- readSegments()
    
    segments$Value <- segments$Value-median(Log2$Value) 
    Log2$Value <- Log2$Value-median(Log2$Value)
    
    cat(' ..processing.\n')
    
    load('allRegions.Rdata')							## These were prepared in TAPS_plot
    #allRegions <- makeRegions(Log2, alf, segments)
    
    ## TODO: Figure out why the following line is trying to load a file that doesn't exist. DCW
    sampleInfo <- load.txt('sampleInfo_TAPS.txt')		## The basis for TAPS copy number analysis
    
    t <- AROMA_findCNs(Log2,alf,allRegions,dmin=0.9,maxCn=maxCn,ceiling=1,shift=sampleInfo$cn2,delta=sampleInfo$delta) ## estimates the Log-R and Allelic Imbalance Ration of all variants up to maxCn
    
    u <- AROMA_setCNs(allRegions,t$int,t$ai,maxCn)			## Assigns copy number variant for all segments
    allRegions$regions <- u$regions
    save.txt(u$merged,file=paste(name,'_segmentCN_merged.txt',sep='')) ## adjacent segments with idendical copy number are merged (except over centromere) and all are saved to a text file
    
    AROMA_scatterPlot_sample(u$regions,t$int,t$ai,xlim=xlim,ylim=ylim,name=name,minseg=minseg,cn=T) ## This plot is saved, marked with TAPS interpretaion of copy number variants. Used for visual quality check.
    for (chr in 1:23) AROMA_scatterChromPlot_Sample(allRegions,Log2,alf,chr,xlim,ylim,name,minseg=5,maxCn=maxCn,width=700,height=1200) # Chromosome-wise plot of all segments with their copy number calls. Used for a vusual quality check.
    cat('..done\n')
    setwd('..')
  }
}
###
AROMA_TAPS_short <- function(directory=getwd(),xlim=c(-1,1),ylim=c(0,1)) {
  ## This function is currently not in use.  	
  setwd(directory)
  subs <- getSubdirs()
  if (is.null(subs)) {
    subs=thisSubdir()
    setwd('..')
  }
  for (i in 1:length(subs)) {
    setwd(subs[i])
    name <- subs[i]
    cat(' ..loading', subs[i])
    Log2 <- AROMA_readLog2()
    alf <- readAlf(localDir)
    segments <- readSegments()
    cat(' ..processing.\n')
    if (is.null(segments)) {
      segments <- AROMA_segment_DNAcopy(Log2)
      save.txt(segments,'segments.txt')
    }
    regs <- AROMA_regsFromSegs(Log2,alf,segments,bin=500,min=1)
    for (chr in 1:24) {
      AROMA_chromPlot_short(chr, Log2, alf, regs, name)
    }
    cat('..done\n')
    setwd('..')
  }
}
###
AROMA_scatterPlot_sample <- function(regions,int=NULL,ai=NULL,xlim,ylim,name,minseg,cn=F) {
  ## TAPS scatter plot of a full sample, used for visual quality control. 
  library(ggplot2)
  regions <- regions[!is.na(regions$imba),]		## Avoids segments that could not be given an Allelic Imbalance Ratio (normally due to having no or almost no SNPs)
  regions <- regions[regions$lengthMB>minseg,]	## May avoid short segments for clarity
  if (cn) {											## If requested, plots the copy numbers (assigned to regions) into the scatter plot
    variants <- unique(regions$fullCN)				## The variants present in this data
    variants_data <- NULL	
    for (v in variants) {							## Extracts the Log-R and Allelic Imbalance Ratio of variants
      t_cn <- strsplit(v,'m')[[1]][1]
      t_int <- as.numeric(int[t_cn][[1]])
      t_ai <- as.numeric(ai[v][[1]])
      variants_data <- rbind(variants_data,c(t_int,t_ai))
    }
    
    png(paste(name,'__TAPS_summary.png',sep='')) ## Plots the whole-sample scatter plot, with copy numbers marked.
    g <- 
      ggplot() + xlab('Log Ratio') + ylab('Allelic Imbalance Ratio') + opts(legend.position="none") + xlim(xlim) + ylim(ylim) +
      geom_point(data=data.frame('x'=regions$log2,'y'=regions$imba,'z'=deChrom_ucsc(regions$Chromosome),'s'=regions$probes),aes(x=x,y=y,size=log10(s), alpha=0.3, colour=ifelse(z==23,'red','black'))) +     
      geom_text(data=data.frame('x'=variants_data[,1],'y'=variants_data[,2],'z'=variants),aes(x=x,y=y,label=z,colour='blue',size=4,alpha=0.7)) +
      scale_colour_identity()
    print(g)
    dev.off()
  } else {
    
    png(paste(name,'__TAPS_summary.png',sep='')) ## If copy number estimates are not provided, plots without copy numbers.
    g <- 
      ggplot() + xlab('Log Ratio') + ylab('Allelic Imbalance Ratio') + opts(legend.position="none") + xlim(xlim) + ylim(ylim) +
      geom_point(data=data.frame('x'=regions$log2,'y'=regions$imba,'z'=deChrom_ucsc(regions$Chromosome),'s'=regions$probes),aes(x=x,y=y,size=log10(s), alpha=0.3, colour=ifelse(z==23,'red','black'))) +        
      scale_colour_identity()
    print(g)
    dev.off()
  }
}
###
AROMA_scatterChromPlot_Sample <- function (allRegions,Log2,alf,chr,xlim,ylim,name,minseg=3,maxCn=12,width=800,height=1000) {
  ## This function builds a TAPS result image of one chromosome. Containing scatter plot, copy number estimates, Log-Ratio and Allele frequency-
  regions <- allRegions$regions						## (allRegions also contains pointers to Log-R and allele frequency of segments)
  regions$imba[regions$lengthMB<minseg] <- NA		## If requested, excludes some very short regions from being plotted
  library(ggplot2)
  qix <- Log2$Chromosome==chrom_ucsc(chr)			## Extracts Log-R, allele freq and segments of this chromosome (chr)
  rix <- alf$Chromosome==chrom_ucsc(chr)
  ix <- regions$Chromosome==chrom_ucsc(chr)
  png(paste(name,'_TAPS_',chrom_ucsc(chr),'.png',sep=''),width=800,height=1000)
  
  siz <- regions$probes								## segments are plotted with respect to their length
  
  p <-												## This is the scatter plot part of the image.
    ggplot() + ylim(ylim) + xlim(xlim) + xlab('Log2(intensity)') + ylab('Allelic imbalance') + opts(legend.position="none") +  
    geom_point(data=data.frame('x'=regions$log2[!ix],'y'=regions$imba[!ix],'col'=sum(ix),'size'=siz[!ix]),aes(x=x,y=y,size=log10(size), alpha=0.5,colour=(-col[1]))) +	## segments on other chromosomes, in grey
    geom_point(data=data.frame('x'=regions$log2[ix],'y'=regions$imba[ix],'col'=1:sum(ix),'size'=siz[ix]),aes(x=x,y=y,size=log10(size), alpha=0.5,colour=col)) +			## segments on this chromosome, color coded on position
    scale_colour_gradient2(low='gray', mid='blue', high='red') + 
    opts(panel.background=theme_rect(colour='white')) 
  
  
  plot_range <- c(min(regions$Start[ix]),max(regions$End[ix]))
  cn <-												## The copy number part of the image
    ggplot()+xlab(NULL)+ylab(NULL) + xlim(plot_range)+ ylim(c(-1,maxCn)) +
    geom_segment(data=data.frame('x'=regions$Start[ix],'y'=regions$mCn[ix],'xend'=regions$End[ix],'yend'=regions$mCn[ix],'col'= -sum(ix)),aes(x=x,y=y,xend=xend,yend=yend, colour=col,size=2,alpha=1)) + ## minor copy numbers are plotted in grey
    geom_segment(data=data.frame('x'=regions$Start[ix],'y'=regions$Cn[ix],'xend'=regions$End[ix],'yend'=regions$Cn[ix],'col'= 1:sum(ix)),aes(x=x,y=y,xend=xend,yend=yend, colour=col,size=2,alpha=1)) +  ## total copy numbers are color coded according to position
    scale_colour_gradient2(low='gray',mid='blue',high='red') + opts(legend.position="none") + opts(plot.margin=unit(c(0,0,0,0),"lines")) +
    opts(panel.background=theme_rect(colour='white'))
  
  col_range <- NULL ## The colors of the Log-R must match segements in other parts of the image
  n <- nrow(regions[ix,])
  for (i in 1:n) {				## for each region...
    region_ix <- which(ix)[i] ## Region's original list position
    region_ix <- allRegions$regionIx$Log2[[region_ix]] ## This region's probe positions in the Log2 vector
    col_range[region_ix] <- i	## put the correct color
  }
  col_range <- col_range[which(qix)]
  
  q <-
    ggplot()+xlab(NULL)+ylab(NULL) + xlim(plot_range)+ ylim(c(-2,3)) +	## The Log-R part of the image
    geom_point(data=data.frame('x'=Log2$Start[qix],'y'=Log2$Value[qix],'c'=col_range),aes(x=x,y=y,colour=c,size=1,alpha=0.1,)) +	## plots all probes
    geom_segment(data=data.frame('x'=regions$Start[ix],'y'=regions$log2[ix],'xend'=regions$End[ix],'yend'=regions$log2[ix],'col'= -max(col_range,na.rm=T)),aes(x=x,y=y,xend=xend,yend=yend, colour=col,size=1,alpha=1)) + ## marks what is a segment
    scale_colour_gradient2(low='black', mid='blue', high='red') + opts(legend.position="none") + opts(plot.margin=unit(c(0,0,0,0),"lines")) +
    opts(panel.background=theme_rect(colour='white'))
  
  col_range <- NULL ## The colors of the alele frequency must match segements in other parts of the image
  n <- nrow(regions[ix,])
  for (i in 1:n) {
    region_ix <- which(ix)[i] # region's original list position
    region_ix <- allRegions$regionIx$alf[[region_ix]] #this region's snp positions in alf vector
    col_range[region_ix] <- i
  }
  col_range <- col_range[which(rix)]
  r <- 
    ggplot()+xlab(NULL)+ylab(NULL)+xlim(plot_range)+	## The allele frequency part of the image
    geom_point(data=data.frame('x'=alf$Start[rix],'y'=alf$Value[rix],'c'=col_range),aes(x=x,y=y,colour=c,size=0.9,alpha=0.3,)) +
    #geom_segment(data=data.frame('chr'=regs$chr[ix],'x'=regs$start[ix],'y'=regs$het[ix],'xend'=regs$end[ix],'col'= -max(col_range)),aes(x=rep(x,2),y=c(0.5+y,0.5-y),xend=rep(xend,2),yend=c(0.5+y,0.5-y), colour=col,size=1,alpha=1)) +
    scale_colour_gradient2(low='black', mid='blue', high='red') + opts(legend.position="none") + opts(plot.margin=unit(c(0,0,0,0),"lines")) +
    opts(panel.background=theme_rect(colour='white'))
  
  
  vplayout(1,131)
  print(p, vp=subplot(1,1:68))	## something very looely based on the Golden Ratio :)
  print(cn, vp=subplot(1,69:89))    
  print(q, vp=subplot(1,90:110))
  print(r, vp=subplot(1,111:131))
  dev.off()
}
###
AROMA_chromPlot_short <- function(chr, Log2, alf, regs, xlim=c(-1.2,1.2),ylim=0:1,name='noname') {
  ## This function is very similar to the one above, except it does not plot copy numbers (they arent known yet), and it uses short
  ## segments that are better suited to a) figure out the sample-specific relationship between Log-R, Allelic Imbalance Ratio and 
  ## allele-specific copy numbers and b) see if segmentation has failed for copy number changes that do not affect total copy number.
  library(ggplot2)
  qix=Log2$Chromosome==chrom_ucsc(chr)
  rix=alf$Chromosome==chrom_ucsc(chr)
  ix=regs$chr==chrom_ucsc(chr)
  if (name=='') name=thisSubdir() 
  png(paste(name,'_',chrom_ucsc(chr),'.png',sep=''),width=800,height=800)
  
  
  p <- 
    ggplot() + ylim(ylim) + xlim(xlim) + xlab('Log2(intensity)') + ylab('Allelic imbalance') + opts(legend.position="none") +  
    geom_point(data=data.frame('x'=regs$logs[!ix],'y'=regs$scores[!ix],'col'=sum(ix)),aes(x=x,y=y,size=1, alpha=0.9,colour=(-col[1]))) + 
    geom_point(data=data.frame('x'=regs$logs[ix],'y'=regs$scores[ix],'col'=1:sum(ix)),aes(x=x,y=y,size=1, alpha=0.9,colour=col)) + 
    scale_colour_gradient2(low='gray', mid='blue', high='red') + 
    opts(panel.background=theme_rect(colour='white')) 
  
  
  plot_range=c(min(Log2$Start[qix]),max(Log2$Start[qix]))
  col_range=regs$key1[qix]-min(regs$key1[qix],na.rm=T)
  q <-
    ggplot()+xlab(NULL)+ylab(NULL) + xlim(plot_range)+ ylim(c(-2,3)) +
    geom_point(data=data.frame('x'=Log2$Start[qix],'y'=Log2$Value[qix],'c'=col_range),aes(x=x,y=y,colour=c,size=1,alpha=0.1,)) +
    geom_segment(data=data.frame('chr'=regs$chr[ix],'x'=regs$start[ix],'y'=regs$logs[ix],'xend'=regs$end[ix],'yend'=regs$logs[ix],'col'= -max(col_range,na.rm=T)),aes(x=x,y=y,xend=xend,yend=yend, colour=col,size=1,alpha=1)) +
    scale_colour_gradient2(low='black', mid='blue', high='red') + opts(legend.position="none") + opts(plot.margin=unit(c(0,0,0,0),"lines")) +
    opts(panel.background=theme_rect(colour='white'))
  
  col_range=regs$key2[rix]-min(regs$key2[rix],na.rm=T)
  r <- 
    ggplot()+xlab(NULL)+ylab(NULL)+xlim(plot_range)+
    geom_point(data=data.frame('x'=alf$Start[rix],'y'=alf$Value[rix],'c'=col_range),aes(x=x,y=y,colour=c,size=0.9,alpha=0.3,)) +
    #geom_segment(data=data.frame('chr'=regs$chr[ix],'x'=regs$start[ix],'y'=regs$het[ix],'xend'=regs$end[ix],'col'= -max(col_range)),aes(x=rep(x,2),y=c(0.5+y,0.5-y),xend=rep(xend,2),yend=c(0.5+y,0.5-y), colour=col,size=1,alpha=1)) +
    scale_colour_gradient(low='blue', high='red') + opts(legend.position="none") + opts(plot.margin=unit(c(0,0,0,0),"lines")) +
    opts(panel.background=theme_rect(colour='white'))
  
  
  vplayout(1,110)
  print(p, vp=subplot(1,1:68))
  print(q, vp=subplot(1,69:89))    
  print(r, vp=subplot(1,90:110))
  dev.off()
}
###
AROMA_regsFromSegs <- function (Log2,alf, segments, bin=200,min=1) {
  ## This function builds short segments and calcualtes their average Log-R and Allelic Imbalance Ratio.
  rownames(Log2)=1:nrows(Log2)
  rownames(alf)=1:nrows(alf)
  regs=list('chr'=NULL,'start'=NULL,'end'=NULL,'logs'=NULL,'scores'=NULL,'het'=NULL,'hom'=NULL,'probes'=NULL,'snps'=NULL,'key1'=rep(NA,nrow(Log2)),'key2'=rep(NA,nrow(alf)))
  n=nrow(segments)
  s_check=NULL
  for (c in 1:n) { ## for every segment
    tlog=Log2[Log2$Chromosome==segments$Chromosome[c],] ## Log-R on this chromosome
    talf=alf[alf$Chromosome==segments$Chromosome[c],] ## Allele Freq on this chromosome
    tlog=tlog[(tlog$Start>=segments$Start[c])&(tlog$Start<segments$End[c]),] ## Log-R on this segment
    talf=talf[(talf$Start>=segments$Start[c])&(talf$Start<segments$End[c]),] ## Allele Freq on this segment
    
    tsnps=nrow(talf)	## number of snps and probes in this segment
    tprobes=nrow(tlog)
    tnregs=max(trunc(tsnps/bin),1) ## Further split into how many (a minimum of 1) subsegments
    tcuts=segments$Start[c] ## The first start pos
    tlength=segments$End[c]-segments$Start[c]	## Length of this whole segment
    for (i in 1:tnregs) tcuts = c(tcuts, tcuts[1]+round(i/tnregs*tlength)) ## Break the segment at these positions
    
    for (r in 1:(tnregs)) {	## build the subsegments
      regs$chr=c(regs$chr,as.character(segments$Chromosome[c]))	## Chromosome
      s_=tcuts[r]													## Start
      e_=tcuts[r+1]												## End
      thisalf=talf[(talf$Start>=s_)&(talf$Start<=e_),]			## get the Log-R values
      thislog=tlog[(tlog$Start>=s_)&(tlog$Start<=e_),]			## and the allele frequency
      regs$key1[as.integer(rownames(thislog))]=length(regs$chr)	## store their positions for fast access during plotting
      regs$key2[as.integer(rownames(thisalf))]=length(regs$chr)	## --"--
      regs$logs=c( regs$logs, mean(thislog$Value) )				## store average log ratio of this segment
      regs$probes=c(regs$probes,nrow(thislog))					## store number of probes
      regs$snps=c(regs$snps,nrow(thisalf))						## store number of bi-allelic probes (SNPs)
      #regs$or_seg=c(regs$or_seg,c)	
      regs$start=c(regs$start,s_)									## store start and end positions
      regs$end=c(regs$end,e_)
      
      if (nrow(thisalf)>min) {									## Time to calculate Allelic Imbalance Ratio (if enough SNPs)
        t1=sort( abs(thisalf$Value-0.5) )							## distance from middle (het) in the allele freq pattern, ascending
        if (length(unique(t1))>3) {								## do not attempt clustering with too few snps
          xx=NULL
          try(xx <- kmeans(t1, 2),silent=T)							## Attempt k-means (Hartigan-Wong: has proven very stable)
          if (!is.null(xx)) if (min(xx$size) > 0.05*max(xx$size)) {	## On some occations data quality is poor, requiring 5%+ heterozygous SNPs avoids most such cases.
            xx=xx$centers
          } else xx=NA
        } else xx=NA      
      } else xx=NA
      if (is.na(xx)) xx=0:1
      regs$scores=c(regs$scores, min(xx)/max(xx) )				## Allelic Imbalance Ratio = inner / outer cluster.
      regs$het=c(regs$het, min(xx))								## $het and $hom are no longer in use.
      regs$hom=c(regs$hom, max(xx))
    }
  }
  return (regs)
}
###
AROMA_segment_DNAcopy <- function(Log2) {
  ## If segmentation is required, DNAcopy is a good choice. Must be installed. 
  library(DNAcopy)
  cat('Using DNAcopy to create segments:\n')
  segs=NULL
  chroms=c(as.character(1:22),'X','Y')
  chroms=paste('chr',chroms,sep='')
  for (c in 1:24) { # segment chromosome
    tlog=Log2[Log2$Chromosome==chroms[c],]	## Log-R of this chromosome
    if (nrow(tlog)>0) {						## (ChrY may be absent)
      cnaObject=segment(smooth.CNA(CNA(tlog$Value, rep(c,nrow(tlog)), tlog$Start, data.type='logratio',sampleid=paste('chr',c))), undo.splits='sdundo', undo.SD=1) ## Segments this chromosome
      segs=rbind(segs,cnaObject$output)		## Add result to data frame
    }
  }
  colnames(segs)=c('x','Chromosome','Start','End','Markers','Value')
  segs$Chromosome=chrom_ucsc(segs$Chromosome)	## Adjust chromosome names to standard
  #save.txt(segs[,2:6],file='_segments.txt')
  return(segs[,2:6]) 
}
###
AROMA_readLog2 <- function() {
  ## This function reads Log-ratio from the file "probes.txt" which must be present in the current directory.
  Log2=NULL
  try( Log2 <- read.csv(file='probes.txt',header=T,sep='\t'), silent=T)
  if (!is.null(Log2)) cat(' ..found probes.txt')
  if (is.null(Log2)) {
    try( Log2 <- read.csv(file='_probes.txt',header=T,sep='\t'), silent=T)
    if (!is.null(Log2)) cat(' ..found _probes.txt')
  }
  ## This code was used if Log-R must be read from .CNCHP file (Affymetrix Genotyping Console or APT). NOT currently supported downstream as .CNCHP is lacks allele-specific information for Affy 250k/500k
  if (is.null(Log2)) {   
    dir=dir()
    ix=grep('cnchp',dir,ignore.case=TRUE)
    if (length(ix)==1)
      cat(' ..found',dir[ix])
    library(affxparser)
    temp=readCcg(dir[ix])
    temp=temp$dataGroups$MultiData$dataSets$CopyNumber$table
    temp$Chromosome[temp$Chromosome==24]=23 # strange numbers in cnchp
    temp$Chromosome[temp$Chromosome==25]=24
    Log2 = data.frame('Chromosome'=chrom_ucsc(temp$Chromosome), 'Start'=temp$Position, 'End'=temp$Position, 'Value'=temp$Log2Ratio)
    save.txt(Log2, file='_probes.txt')
    cat(' ..wrote _probes.txt')
    adif=NULL
    try(
      adif <- temp$Allele,
      silent=T
    )      
    if (!is.null(adif)) {
      ix= !is.nan(adif)
      alf=data.frame('Chromosome'=chrom_ucsc(temp$Chromosome[ix]), 'Start'=temp$Position[ix], 'End'=temp$Position[ix], 'Value'=0.5+adif[ix])
      save.txt(alf, file='_snps.txt')
      cat(' ..wrote _snps.txt')
    }
  }
  return (Log2) 
}
###
readAlf <- function(localDir=NULL) {
  ## This funciton reads allele frequency [B/(A+B)] from the file 'snps.txt', which must be present in the current directory.
  alf=NULL
  try( alf <- read.csv(file='snps.txt',header=T,sep='\t'), silent=T)
  if (!is.null(alf)) cat(' ..found snps.txt')
  if (is.null(alf)) {
    try( alf <- read.csv(file='_snps.txt',header=T,sep='\t'), silent=T)
    if (!is.null(alf)) cat(' ..found _snps.txt')
  }
  return (alf) 
}
###
readSegments <- function() {
  ## This function reads segment information in 'segments.txt' which must be present in the current folder.
  ## The author recommends SNP-rank segmentation (NEXUS) or another CBS such as that in DNACopy.
  ## Using a HMM is not recommended unless you have a homogenous, diploid sample. (And then there is more user-friendly software anyway.)
  segments=NULL
  try( segments <- read.csv(file='segments.txt',header=T,sep='\t'),silent=T)
  if (!is.null(segments)) {
    cat(' ..found segments.txt')
  }
  if (is.null(segments)) { 
    try( segments <- read.csv(file='_segments.txt',header=T,sep='\t'),silent=T)
    if (!is.null(segments)) {
      cat(' ..found _segments.txt')
    }
  }
  return (segments)
}
###
makeRegions <- function(Log2, alf, segments,dataType='Nexus') {
  ## makeRegions is similar to "regsfromsegs" except regions are not subdivided before calculation of mean Log-R and Allelic Imbalance Ratio.
  regions=segments
  regions$Chromosome=as.character(segments$Chromosome)			## Chromosome
  regions$lengthMB=round((regions$End-regions$Start)/1000000,3)	## length in megabases
  regions$probes=0												## # of probes
  regions$snps=0													## # of bi-allelic probes (SNPs)
  #regions$log2=round(regions$Value,4)
  regions$imba=NA													## Allelic Imbalance Ratio
  regionIx=NULL													## Not currently used
  for (i in 1:nrows(regions)) {
    log2temp=which(equals(Log2$Chromosome,regions$Chromosome[i])) ## index of Log-R (current chrom)
    alftemp=which(equals(alf$Chromosome,regions$Chromosome[i]))	## index of Allele frequency (current chrom)
    
    log2temp=log2temp [Log2$Start[log2temp]>=regions$Start[i] & Log2$Start[log2temp]<regions$End[i]]	## index of Log-R (current segment)
    alftemp=alftemp [alf$Start[alftemp]>=regions$Start[i] & alf$Start[alftemp]<regions$End[i]]		## index of Allele frequency (current segment)
    regions$probes[i]=length(log2temp)
    regions$snps[i]=length(alftemp)
    regionIx$Log2[[i]]=log2temp			## indexes of Log-R and Allele Frequency are saved for future use (plot color coding)
    regionIx$alf[[i]]=alftemp
    log2temp=Log2$Value[log2temp]
    regions$log2[i]=median(log2temp)
    alftemp=alf$Value[alftemp]
    if (length(alftemp)>3) {				## Prepare to calculate Allelic Imbalance Ratio
      if (dataType=="Nexus") t1=sort( abs(alftemp-0.5) ) ## distance from middle (het), ascending
      if (dataType=="CNCHP") t1=sort( abs(alftemp)-0 ) ## This is not currently in use, was intended for analysis of SNP6 .CNCHP data.
      if (length(unique(t1))>5) { ## Avoid calculating Allelic Imbalance Ratio unless there are several different values to cluster
        xx=kmeans(t1, 2)			## This part is nearly identical to that of 'regsfromsegs()'
        if (min(xx$size) > 0.05*max(xx$size)) {
          xx=xx$centers
        } else xx=NA	
      } else xx=NA
      regions$imba[i]=round( min(xx)/max(xx) ,2)
    }
  }
  return(list('regions'=regions,'regionIx'=regionIx))
}
###
equals <- function(a,b) {
  ## a helper function in case factors are compared
  a <- as.character(a)
  b <- as.character(b)
  return (a==b)
}
###
nrows <- function(data) dim(data)[1] ## don't ask me about this one.
###
load.txt <- function(file) { ## just to be rid of "sep"
  read.csv(file,sep='\t') 
}
###
save.txt <- function(data,file='temp', row.names=F, col.names=T, append=F) {   ## simplified for convenience
  write.table(data,file=file,sep='\t',quote=F,row.names=row.names, col.names=col.names, append=append)
}
###
deChrom_ucsc <- function(data) {
  ##  chroms in 1:24 and chr1:chrY format conversion
  keep_index=NULL
  CHR=rep(0,length(data))
  existing=unique(data)
  for (i in 1:length(existing))   {
    temp=strsplit(as.character(existing[i]),'_')
    temp=strsplit(temp[[1]],'chr')[[1]][2]
    if (temp=='X') temp='23'
    if (temp=='_X') temp='23'
    if (temp=='Y') temp='24'
    CHR[data==existing[i]]=as.integer(temp)
  }
  return (CHR)
}
###
chrom_ucsc <- function(data) {
  ## same as above, backwards
  n=length(data)
  out=rep("",n)
  for (i in 1:n) {
    temp=as.character(data[i])
    if (temp=='23') temp='X'
    if (temp=='24') temp='Y'
    out[i]=paste('chr',temp, sep='')
  }
  return(out)
}
###
subplot <- function(x, y) viewport(layout.pos.col=x, layout.pos.row=y) ## for subplots
vplayout <- function(x, y) {		
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(y,x)))
}
###
getSubdirs <- function() {
  ## Reach for the subdirectories of current directory.
  dir=dir()
  here=getwd()
  n=length(dir)
  subdir=F
  dir_ix=rep(FALSE,n)
  for (i in 1:n) {
    try(setwd(dir[i]), silent=T)
    if (getwd()!=here) { # in another directory
      setwd('..')
      subdir=T
      dir_ix[i]=T
    }
  }
  if (subdir) return (dir[dir_ix])
  return (NULL)
}
###
thisSubdir <- function() { 
  ## get the name of this subdirectory
  here=getwd()
  here=strsplit(here,'/')[[1]]
  here=here[length(here)]
  return(here)
}
###
mySorter <- function(data) {
  ## sort data on $Chromosome and $Start
  temp=data[order(data$Start),]
  t=deChrom_ucsc(temp$Chromosome)
  temp=temp[!is.na(t),]
  t=t[!is.na(t)]
  temp=temp[order(t),]
  return (temp)
}
###
# probesInRegions <- function(regions,probes) {
#   ## This function is not in use.
#   rchr=as.character(regions$Chromosome)
#   pchr=as.character(probes$Chr)
#   nprobes=0
#   
#   for (i in 1:length(rchr)) {
#     t=strsplit(rchr[i],'chr')[[1]][2]
#     here=probes[(t==pchr),] #finding probe subset on this chr
#     matching=(here$Position>=regions$Start[i]) & (here$Position<=regions$End[i])
#     nprobes[i]=sum(matching)
#   }
#   return(nprobes)
# }
###
getOverlap <- function(s1,e1,s2,e2) { ## This is not currently used anywhere - it is destined for group studies and pileups.
  return(max((min(e1,e2)-max(s1,s2)),0))
}
###
weightedMedian <- function(data,weights) {
  try( if (length(weights)==0) return (NULL), silent=T)
  try( return (median(rep(data,weights))), silent=T)
  return (NULL)
}
weightedMean <- function(data,weights) {
  if (length(weights)==0) return (NULL)
  return (mean(rep(data,weights)))
}
###
# myMeans <- function(data,k) { ## Not currently used
#   best <- 0
#   clus <- list()
#   for (i in 1:10) {
#     clus[[1]] <- kmeans(data,k)
#     var <- var(clus[[1]]$centers)
#     if (var>best) {
#       ix <- i
#       best <- var
#     }
#   }
#   return (clus[[ix]])
# }
is.even <- function(data) {
  return (trunc(data/2)*2 == data)
}
is.odd <- function(data) {
  return (trunc(data/2)*2 != data)
}
neighbour <- function(data) {
  n <- length(data)
  dif <- data[2:n]-data[1:(n-1)]
  return (mean(dif))
}
is.autosome <- function(vector) {
  temp <- deChrom_ucsc(vector)
  return (temp<23)
}

## 
AROMA_getEstimates <- function(Log2,alf,regions) {
  ## This function is called to figure out copy number two and the distance to copy number 1 or 3. Curently rather stable unless 
  ## a) normal cn2 (copy-number-two) is missing or b) sample is tetraploid+ or c) Normal cell content is very high.
  ## Therefore inspection/correction of estimates is NECESSARY before proceeding with TAPS_call. 
  ## This will be improved further.
  median <- median_low <- median_up <- NULL
  temp <- regions[(is.autosome(regions$Chromosome)&regions$lengthMB>3)&(!is.na(regions$imba)),]	## Use only long regions >3MB here, for stability
  median <- weightedMedian(temp$log2,temp$probes)				## In many samples, this is copy number 2
  ix <- temp$log2<(median-0.1)
  median_low <- ifelse (sum(ix)>0, weightedMedian(temp$log2[ix],temp$probes[ix]), median-0.4)	## Which implies this is copy number 1
  ix <- (temp$log2 > (median+0.1) ) & (temp$log2 < median + (median-median_low))
  median_up <- ifelse (sum(ix)>0, weightedMedian(temp$log2[ix],temp$probes[ix]), median+0.3) ## And this is copy number 3
  
  t1 <- temp$imba[abs(temp$log2-median_low)<0.1]
  t2 <- temp$imba[abs(temp$log2-median)<0.1]
  t3 <- temp$imba[abs(temp$log2-median_up)<0.1]
  ## if lowest AIR (median) < median_up < median_low we are ok.	    	    
  if ((min(t1)>min(t3))&(min(t3)>min(t2))) return(data.frame('cn2'=median,'delta'=median-median_low))
  
  ## "median" might be copy number 3. Then lowest AIR (median_low) < median > median_up
  if ((min(t1)<min(t2))&(min(t3)<min(t2))) return(data.frame('cn2'=median_low,'delta'=1.2*(median-median_low)))
  
  ## "median" might occationally be copy number 4. Then lowest AIR (median_low) > median < median_up but cn2 is still likely median_low.
  if ((min(t1)>min(t2))&(min(t3)>min(t2))) return(data.frame('cn2'=median_low,'delta'=(median-median_low)))
  
  ## "median" might be copy number 3 with the special case that (median_low) > median > median_up
  if ((min(t1)>min(t2))&(min(t3)<min(t2))) return(data.frame('cn2'=median_low,'delta'=1.2*(median-median_low)))
  
  return (data.frame('cn2'=median,'delta'=median-median_low)) ## default.
}
###
AROMA_findCNs <- function(Log2,alf,allRegions,name=thisSubdir(),dmin=0.9,maxCn=10,ceiling=1,shift=NULL,delta=NULL) {
  ## This function takes an estimate of the Log-R of copy number two (shift) and the difference in log-R between copy numbers 2 and 1 (delta)
  ## (3 and 2 works too). Then, the Log-R and Allelic Imbalance Ratio of all possible copy number variants up to maxCn are estimated from
  ## the Log-R and Allelic Imbalance Ratio of all the segments. This function will NOT be useful unless there is already a solid estimate 
  ## of 'shift' and 'delta'. See the previous function.
  tix=NULL	#temporary index
  int=NULL	## contains Log-R estimate of each (total) copy number
  ai=NULL		## contains Allelic Imbalance Ratio estimate of each copy number variant.
  regions <- allRegions$regions
  regions <- regions[(is.autosome(regions$Chromosome)&regions$lengthMB>1)&(!is.na(regions$imba)),] ## will use these regions
  
  ## likely cn2 regions sit within delta/3 from shift.
  expectedAt <- shift
  tix$cn2 <- abs(regions$log2 - expectedAt) < (delta/3)	## index of likely cn2 regions
  temp <- regions[tix$cn2,]								## cn2 regions
  med <- weightedMedian(temp$log2,temp$probes)			## improved value of Log-R at cn2 (returns NULL if theres nothing there)
  int$cn2 <- ifelse(!is.null(med),med,expectedAt)			## saved to int.
  
  ## likely cn1 regions sit at about cn2 - delta:
  d <- delta												## the (Log-R) distance to cn1
  expectedAt <- int$cn2-d									## cn1 is expected here
  tix$cn1 <- abs(regions$log2 - expectedAt) < (d/3)		## index of likely cn1 regions
  temp <- regions[tix$cn1,]
  med <- weightedMedian(temp$log2,temp$probes)
  int$cn1 <- ifelse(!is.null(med),med,expectedAt)
  
  ## likely cn0 regions sit below cn1 - delta:
  d <- int$cn2-int$cn1
  expectedAt <- int$cn1-d									
  tix$cn0 <- regions$log2 < expectedAt
  temp <- regions[tix$cn0,]
  med <- weightedMedian(temp$log2,temp$probes)
  int$cn0 <- ifelse(!is.null(med),med,expectedAt)
  
  ## likely cn3 regions sit at about cn2+delta*dmin (dmin is about 0.7-0.9, the amount by which the distance between consecutive copy numbers diminish)
  d <- delta*dmin
  expectedAt <- int$cn2+d
  tix$cn3 <- abs(regions$log2 - expectedAt) < (d/3)
  temp <- regions[tix$cn3,]
  med <- weightedMedian(temp$log2,temp$probes)
  int$cn3 <- ifelse(!is.null(med),med,expectedAt)
  
  ## cn4 follows at ...
  d <- dmin*(int$cn3-int$cn2) 
  expectedAt <- int$cn3+d
  tix$cn4 <- abs(regions$log2 - expectedAt) < (d/4)
  temp <- regions[tix$cn4,]
  med <- weightedMedian(temp$log2,temp$probes)
  int$cn4 <- ifelse(!is.null(med),med,expectedAt)
  
  ## generalized for higher cns
  for (cn in 5:maxCn) {
    thisCn <- paste('cn',cn,sep='')
    prevCn <- paste('cn',cn-1,sep='')
    pprevCn <- paste('cn',cn-2,sep='')
    d <- dmin*(int[prevCn][[1]]-int[pprevCn][[1]])
    expectedAt <- int[prevCn][[1]]+d
    tix[[thisCn]] <- abs(regions$log2 - expectedAt) < (d/5)
    temp <- regions[tix[thisCn][[1]],]
    med <- weightedMedian(temp$log2,temp$probes)
    int[thisCn] <- ifelse(!is.null(med),med,expectedAt)
  }
  
  
  ## at cn2, find the variant clusters (normal and CNNLOH)
  ix <- (abs(regions$log2 - int$cn2) < 0.2*(int$cn3-int$cn2) ) # taking only closely-matching segments
  data <- regions[ix,]
  data <- data[!is.na(data$imba),] # ...with a calculated allelic imbalance.
  
  attempt <- 0; while (T) {
    attempt <- attempt+1
    centers <- sort(kmeans(rep(data$imba,data$snps),2)$centers)   ## May with improvement work on even lower T%
    if (((centers[1]<0.2)&(centers[2]<0.8))&(centers[2]>0.3)) break
    if (attempt > 20) { # normal or CNNLOH may be missing
      if (weightedMedian(data$imba,data$snps)>0.3) { #normal seems to be missing. set it to 0.12
        centers <- c(0.12, weightedMedian(data$imba,data$snps))
        break
      }else if (weightedMedian(data$imba)<0.25) { #CNNLOH seems to be missing. set it while looking at cn3
        centers <- c(weightedMedian(data$imba,data$snps),NA)
        break
      }else {	## Neither seems to work. 
        centers <- c(0.12,0.56) ## This may suffice for downstream estimates to work out.
        break
      }
    }
  }
  
  ai$cn2m1 <- centers[1]
  ai$cn2m0 <- centers[2]
  
  ## for cn1 (and 0)
  ix <- (abs(regions$log2 - int$cn1) < 0.2*(int$cn2-int$cn1) )
  data <- regions[ix,]
  data <- data[!is.na(data$imba),]
  expectedAt <- (ai$cn2m0+ai$cn2m1)*3/5	## Decent estimate.
  med <- weightedMedian(data$imba,data$snps) ## Average allelic imbalance weighted on snp count
  ai$cn1m0 <- ifelse (!is.null(med),med,expectedAt)	## Will be NA if there was no CNNLOH
  ai$cn0m0 <- NA #unimportant
  
  ## for cn3:
  ix <- (abs(regions$log2 - int$cn3) < 0.2*(int$cn4-int$cn3) )
  data <- regions[ix,]
  data <- data[!is.na(data$imba),]
  
  if (!is.na(ai$cn2m0)) { # CNNLOH is not missing
    range <- ai$cn2m0-ai$cn2m1 #the distance between 2normal and CNNLOH
    # get the 3(1) regions: 
    expectedAt <- ai$cn2m1+range/3 # this is an approx if cn3m1 is absent.
    
    ix <- (data$imba<ai$cn2m0) & (data$imba>ai$cn2m1) # take regions with less AI than cn2m0 but more than cn2m1
    med <- weightedMedian(data$imba[ix],data$snps[ix]) # average allelic imbalance weighted on snp count
    ai$cn3m1 <- ifelse (!is.null(med),med,expectedAt)
    
    # now for cn3m0
    expectedAt <- ai$cn3m1 + range*dmin # approx of cn3m0
    
    ix <- (abs(data$imba-expectedAt) / abs(data$imba-ai$cn3m1)) < 0.5  # take those much closer to exp than to cn3m1
    med <- weightedMedian(data$imba[ix],data$snps[ix]) 
    try (if (med<ai$cn2m0) med <- NULL, silent=T)	## in case of heterogeneities or strange signals, we disallow this effect.
    ai$cn3m0 <- ifelse (!is.null(med),med,expectedAt)
  } else { ## (if CNNLOH was missing, use a different strategy.)
    data <- data[data$imba>ai$cn2m1,] # avoid 2(0) / 4(0) heterogeneity
    data <- data[data$imba<0.85,] # avoid highly-hom (clus-errors)
    #centers <- sort(kmeans(rep(data$imba,data$snps),2)$centers)## bad idea if either is missing! (common on low-mut samples)
    ix <- data[data$imba<0.4,]
    med <- weightedMedian(data$imba[ix],data$snps[ix])
    ai$cn3m1 <- ifelse (!is.null(med),med,ai$cn2m1+0.1)
    ix <- data[data$imba>0.4,]
    med <- weightedMedian(data$imba[ix],data$snps[ix])
    ai$cn3m0 <- ifelse (!is.null(med),med,ai$cn2m1+0.35)
    ai$cn2m0 <- ai$cn2m1+(ai$cn3m0-ai$cn3m1) # the CNNLOH estimate from CN3
    #print(cat('Warning: cn2m0 missing, estimated cn3m1,cn3m0 centers at',ai$cn3m1,ai$cn3m0,'\n'))
  }
  
  if (is.na(ai$cn1m0)) ai$cn1m0 <- (ai$cn3m1+ai$cn2m0)/2 ## If deletions were missing, place an estimate from cn3 
  
  ## now for cn4
  ix <- abs(regions$log2 - int$cn4) < 0.2*(int$cn4-int$cn3) 
  data <- regions[ix,]
  data <- data[!is.na(data$imba),]
  
  # get the 4(2) regions: they are at AI about cn2m1
  expectedAt <- ai$cn2m1 # this is a good approx 
  
  ix <- data$imba<ai$cn3m1 # let all below cn3m1 in
  med <- weightedMedian(data$imba[ix],data$snps[ix]) 
  ai$cn4m2 <- ifelse (!is.null(med),med,expectedAt)
  
  data <- data[!ix,] # coutinue with the remaining
  # now 4(1) has less ai than 3(0) -> sits at about cn3m1+(cn3m1-cn2m1).
  expectedAt <- ai$cn3m1+(ai$cn3m1-ai$cn2m1)
  ix <- abs(data$imba-expectedAt)<abs(data$imba-ai$cn3m0) # take those closer to exp than to 3(0)
  med <- weightedMedian(data$imba[ix],data$snps[ix]) 
  ai$cn4m1 <- ifelse (!is.null(med),med,expectedAt)
  
  # now 4(0) is at about 3(0) + [0.9 - 3(0)] * [3(0)-2(0)] / [0.9-2(0)] 
  expectedAt <- ai$cn3m0 + (ceiling-ai$cn3m0)*(ai$cn3m0-ai$cn2m0)/(ceiling-ai$cn2m0)
  ix <- abs(data$imba-expectedAt) < 0.2*(expectedAt-ai$cn4m2) # we take those very close to exp
  med <- weightedMedian(data$imba[ix],data$snps[ix])
  try (if (med<ai$cn3m0) med <- NULL, silent=T)	## in case of heterogeneities or strange signals, we disallow this effect.
  ai$cn4m0 <- ifelse (!is.null(med),med,expectedAt)
  
  ## generalization for higher copy numbers. it gets complicated.
  for (cn in 5:maxCn) {
    thisCn <- paste('cn',cn,sep='')
    prevCn <- paste('cn',cn-1,sep='')
    pprevCn <- paste('cn',cn-2,sep='')
    
    ix <- (abs(regions$log2 - int[thisCn][[1]]) < 0.2*(int[thisCn][[1]]-int[prevCn][[1]]) )
    data <- regions[ix,]
    data <- data[!is.na(data$imba),]
    data <- data[data$lengthMB>3,] # long regions for safety
    
    ## try to find variants, starting with LOH
    # LOH such as 5(0)
    m <- 0
    thisVariant=paste(thisCn,'m',0,sep='')
    c4m0 <- ai[paste(prevCn,'m',m,sep='')][[1]] # relative naming for clarity
    c3m0 <- ai[paste(pprevCn,'m',m,sep='')][[1]]
    
    expectedAt <- ceiling-((ceiling-c4m0)*(ceiling-max(c4m0,c3m0))/(ceiling-min(c3m0,c4m0)))
    #ai[thisVariant] <- expectedAt
    # we then take those close to exp (within delta[3(0),4(0)])
    ix <- abs(data$imba-expectedAt) < abs(c4m0 - c3m0) 
    med <- weightedMedian(data$imba[ix],data$snps[ix]) 
    try (if (med<ai$cn4m0) med <- NULL, silent=T)	## in case of heterogeneities or strange signals, we disallow this effect.
    ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
    
    ## then from balanced to less balanced
    minorVariants=trunc(cn/2):1 
    first <- T
    for (m in minorVariants) {
      thisVariant=paste(thisCn,'m',m,sep='')
      if (m==cn/2) {
        # We have balanced variant, so rather easy.
        expectedAt <- ai[paste(pprevCn,'m',m-1,sep='')][[1]] # this is a good approx, balanced at cn-2
        ix <- data$imba<ai[paste(prevCn,'m',m-1,sep='')][[1]] # let all below (cn-1, mcn-1) in
        med <- weightedMedian(data$imba[ix],data$snps[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
        first <- F
      } else if (first) { 
        # its not balanced but its the most balanced of the unbalanced. something like 5(2)
        expectedAt <- 0.5*( ai[paste(prevCn,'m',m,sep='')][[1]] + ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) # that means between 4(2) and 3(1)
        ix <- abs(data$imba-expectedAt) < ( ai[paste(prevCn,'m',m,sep='')][[1]] - ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) /3 # let all "between" 4(2) and 3(1) in
        med <- weightedMedian(data$imba[ix],data$snps[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
        first <- F
      } else {
        # not the most balanced unbalanced variant, for example 5(1):
        expectedAt <- ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]] + (minorVariants[1]-m) * (ai[paste(thisCn,'m',0,sep='')][[1]] - ai[paste(thisCn,'m',minorVariants[1],sep='')][[1]]) / trunc(cn/2) # 5(2) + (which)* 5(0)-5(2) /(n)
        ix <- abs(data$imba-expectedAt) < ( ai[paste(prevCn,'m',m,sep='')][[1]] - ai[paste(pprevCn,'m',m-1,sep='')][[1]] ) /3 # let all "between" 4(1) and 3(0) in
        med <- weightedMedian(data$imba[ix],data$snps[ix]) 
        ai[thisVariant] <- ifelse (!is.null(med),med,expectedAt)
      } 
    } # done with minor variants
  } # done with copy numbers
  
  return(list('int'=int,'ai'=ai))
}
###
AROMA_setCNs <- function(allRegions,int,ai,maxCn=12) {
  ##  Assign total and minor copy numbers to all segments.
  regions <- allRegions$regions[,-4]	## This time, work on all segments available.
  
  Cn <- NULL			## Total copy number
  mCn <- NULL			## Minor allele copy number
  fullCN <- NULL		## Variant label. ('cnXmY')
  
  intDist <- NULL		## distance to certain Log-R
  imbaDist <- NULL	## distance to certain allelic imbalance
  
  ## give each segment the best matching CN and mCN
  for (i in 1:nrow(regions)) {
    # set total copy number
    distance <- Inf
    for (cn in 0:maxCn) {
      t_int <- int[paste('cn',cn,sep='')][[1]]	## get Log-R of particular cn from 'int'
      t_dis <- abs(regions$log2[i]-t_int)			## distance to that particular cn
      if (t_dis < distance) {						## nearest so far, save.
        distance <- t_dis -> intDist[i]
        Cn[i] <- cn
      }
    }
    
    # Y makes CN0 sit very low, this is a fix on non-Y Cn0
    #if ((Cn[i] == 1)&(intDist[i] > (int$cn2-int$cn1))) Cn[i] <- 0 ## currently not needed.
    
    # set minor CN
    distance <- Inf
    if (Cn[i]<=1) {
      mCn[i] <- 0
    } else if (!is.na(regions$imba[i])) for (m in 0:trunc(Cn[i]/2)) {
      t_ai <- ai[paste('cn',Cn[i],'m',m,sep='')][[1]]
      t_dis <- abs(regions$imba[i]-t_ai)
      if (t_dis < distance) {
        distance <- t_dis -> imbaDist[i]
        mCn[i] <- m
      }
    } else mCn[i] <- NA
    fullCN[i] <- paste('cn',Cn[i],'m',mCn[i],sep='') # Full description
  }
  
  ## label cn1 and cn0 - not currently needed.
  #mCn[Cn==1] <- 0
  #fullCN[Cn==1] <- 'cn1m0'
  #mCn[Cn==0] <- 0
  #fullCN[Cn==0] <- 'cn0m0'
  
  regions$Cn <- Cn
  regions$mCn <- mCn
  regions$fullCN <- fullCN
  
  data <- regions[,c(-7,-8)]
  ## Merge consecutive identical variants
  row <- 1
  while (row < nrow(data)) { # while not on the last row
    #print(row)
    if (((data$Chromosome[row] == data$Chromosome[row+1]) & (data$fullCN[row] == data$fullCN[row+1])) & ( (data$Start[row+1]-data$End[row])<5000 ) ) { ## segments are adjacent with same copy number and not separated by 5kb+ (centromere)
      data$End[row] <- data$End[row+1]
      data$probes[row] <- data$probes[row]+data$probes[row+1]
      data$snps[row] <- data$snps[row]+data$snps[row+1]
      data$lengthMB[row] <- data$lengthMB[row]+data$lengthMB[row+1]
      data <- data[-(row+1),]  # Merges and deletes row
    } else {
      row <- row+1 # or leaves it
    }
  }
  return (list('regions'=regions,'merged'=data))
}


AROMA_CRMAV2_toTAPS <- function(CEL_directory=getwd(), poolRef=TRUE) {
  ## Use AROMA CRMAV2 to generate normalized log-R and allele ratio data. 
  ## AROMA must be installed, setup must be completed correctly, please see
  ## http://www.aroma-project.org/getstarted.
  ## Note that CEL files should reside in the AROMA subfolder rawData/[projectName]/[arrayType].
  ## Run this function from within the folder with the CEL files (or specify the folder as argument).
  ## Sample set average will be used as reference, to run without reference add argument poolRef=FALSE.
  ## Subfolder will be made with TAPS-compatible data. Once this function is finished,
  ## go ahead with AROMA_TAPS_plot().
  
  library("aroma.affymetrix")
  
  setwd(CEL_directory)
  celfiles <- dir()[grep('.CEL',dir())]
  nsamples=length(celfiles)
  samples=NULL
  for (i in 1:nsamples) {
    samples[i]=strsplit(celfiles[i],'.CEL')[[1]][1]
    try(dir.create(samples[i]), silent=T)
  }
  
  chipType=thisSubdir()
  setwd('..')
  projName=thisSubdir()
  setwd('../..')
  aromaDir=getwd()
  
  log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
  
  dataSet <- projName;
  
  options(digits=4)
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
  gi <- getGenomeInformation(cdf)
  si <- getSnpInformation(cdf)
  
  acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
  
  csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf)
  
  acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
  #print(acc)
  csC <- process(acc, verbose=verbose)
  print(csC)
  
  #By using argument target="zero", no reference is required.  Otherwise, the average file will be used as the reference.
  bpn <- BasePositionNormalization(csC, target="zero")
  print(bpn)
  csN <- process(bpn, verbose=verbose)
  print(csN)
  
  
  plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE)
  #plmtot <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE)
  
  if (length(findUnitsTodo(plm)) > 0) {
    units <- fitCnProbes(plm, verbose=verbose)
    units <- fit(plm, verbose=verbose)
  }
  
  #if (length(findUnitsTodo(plmtot)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
  #units <- fitCnProbes(plmtot, verbose=verbose)
  #str(units)
  # Fit remaining units, i.e. SNPs (~5-10min/array)
  #units <- fit(plmtot, verbose=verbose)
  #str(units)
  #}
  
  
  ces <- getChipEffectSet(plm)
  #cestot <- getChipEffectSet(plmtot)
  
  fln <- FragmentLengthNormalization(ces, target="zero")
  #flntot <- FragmentLengthNormalization(cestot, target="zero")
  
  cesN <- process(fln, verbose=verbose)
  #cesNtot <- process(flntot, verbose=verbose)
  print(cesN)
  
  cdf <- getCdf(cesN) 
  gi <- getGenomeInformation(cdf)
  
  ## postprocess 
  
  for (i in 1:nsamples) {
    
    ce <- getFile(cesN, i)  # Array #3
    ceR <- getAverageFile(cesN)
    
    #cetot <- getFile(cesNtot, i)
    Log2 <- alf <- NULL
    for (chromosome in 1:23) {
      units <- getUnitsOnChromosome(gi, chromosome=chromosome)
      pos <- getPositions(gi, units=units)
      
      theta <- extractTheta(ce, units=units)
      baf <- theta[,2] / (theta[,1]+theta[,2])
      tot=theta[,1]+theta[,2]
      thetaR <- extractTheta(ceR, units=units)
      totR=thetaR[,1]+thetaR[,2]
      
      if (poolRef==FALSE) totR=mean(tot)
      
      tot=log(tot/totR)
      #thetatot <- extractTheta(cetot, units=units)
      
      baf=baf[order(pos)]
      log=tot[order(pos)]
      pos <- pos[order(pos)]
      
      ix <- !is.na(log)
      Log2 <- rbind(Log2, data.frame('Chromosome'=chrom_ucsc(chromosome), 'Start'=pos[ix], 'End'=pos[ix], 'Value'=log[ix]))
      ix <- !is.na(baf)
      alf <- rbind(alf, data.frame('Chromosome'=chrom_ucsc(chromosome), 'Start'=pos[ix], 'End'=pos[ix], 'Value'=baf[ix]))
    }
    setwd(CEL_directory)
    setwd(samples[i])
    save.txt(Log2, 'probes.txt')
    save.txt(alf, 'snps.txt')
    setwd(aromaDir)
  }
  setwd(CEL_directory)
  cat(nsamples, 'samples processed and stored for TAPS analysis. Next step is to run AROMA_TAPS_plot().\n')
}


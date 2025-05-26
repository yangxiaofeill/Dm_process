## repeat masker processing
setwd("~/mode-species/fly/RepeatMasker/")

library(data.table)
library(ggplot2)
library(cowplot)
library(IRanges)
library(dplyr)
library(purrr)

rm(list = ls())

getFreq <- function(theRepeatMasker, GenomeSize){
  chrs <- unique(theRepeatMasker$V5)
  freqRes <- data.frame(chr = "", rType = "", size = 0, percentage = 0.0)
  freqRes <- freqRes[-1, ]
  for(oneChr in chrs){
    message("************ ", oneChr, " ************")
    oneChrRM <- theRepeatMasker[which(theRepeatMasker$V5 == oneChr), ]
    theTypes <- unique(oneChrRM$V11)
    chrSize <- GenomeSize$V2[which(GenomeSize$V1 == oneChr)]
    oneChrFreq <- data.frame(chr = "", rType = "", size = 0, percentage = 0.0)
    oneChrFreq <- oneChrFreq[-1, ]
    for(oneType in theTypes){
      oneChrOneTypeRM <- oneChrRM[which(oneChrRM$V11 == oneType), ]
      oneChrOneTypeSize <- sum(oneChrOneTypeRM$V7 - oneChrOneTypeRM$V6)
      oneP <- oneChrOneTypeSize / chrSize
      message("        ", oneChr, ": ", oneType, ": Size = ", oneChrOneTypeSize, ", Percentage: ", oneP)
      oneChrFreq <- rbind(oneChrFreq, data.frame(chr = oneChr, rType = oneType, size = oneChrOneTypeSize, percentage = oneP))
    }
    oneChrFreq <- oneChrFreq[order(oneChrFreq$size, decreasing = T), ]
    
    freqRes <- rbind(freqRes, oneChrFreq)
  }
  return(freqRes)
}

Msize <- fread("../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)
M_RepeatMasker <- fread("./M.np2.manule.fa.out", skip = 3, stringsAsFactors = F, fill = T)
M_repeatFreq <- getFreq(theRepeatMasker = M_RepeatMasker, GenomeSize = Msize)
write.table(M_repeatFreq, file = "./M.np2.manule.fa.out.freq", quote = F, sep = "\t", row.names = F, col.names = F)

Wsize <- fread("../W.np2.manule.fa.fai", header = F, stringsAsFactors = F)
W_RepeatMasker <- fread("./W.np2.manule.fa.out", skip = 3, stringsAsFactors = F, fill = T)
W_repeatFreq <- getFreq(W_RepeatMasker, Wsize)
write.table(W_repeatFreq, file = "./W.np2.manule.fa.out.freq", quote = F, sep = "\t", row.names = F, col.names = F)

getBins <- function(genomeSize, binSize){
  starts <- seq(1, genomeSize, by = binSize)
  ends <- seq(binSize, genomeSize, by = binSize)
  ends <- c(ends, genomeSize)
  
  theBins <- data.frame(start = starts, end = ends)
  return(theBins)
}

getBinPvalues <- function(oneTypeRM, genomeSize, binSize){
  theBins <- getBins(genomeSize, binSize)
  binRange <- IRanges(start = theBins$start, end = theBins$end)
  repeatRange <- IRanges(start = oneTypeRM$V6, end = oneTypeRM$V7)
  
  xmins <- oneTypeRM$V6
  xmaxs <- oneTypeRM$V7
  total_coverage <- sum(xmaxs - xmins + 1)
  
  # calculate the overlap length
  overlaps <- findOverlaps(binRange, repeatRange)
  intersect_lengths <- width(pintersect(
    binRange[queryHits(overlaps)],
    repeatRange[subjectHits(overlaps)]
  ))
  coverage_lengths <- tapply(intersect_lengths, queryHits(overlaps), sum)
  theBins$k <- 0
  theBins$k[as.integer(names(coverage_lengths))] <- coverage_lengths
  
  ## calculate the p-value
  window_length = theBins$end - theBins$start + 1
  coverage_density <- total_coverage / genomeSize
  E <- coverage_density * window_length
  p_value <- pbinom(q = theBins$k, size = window_length, prob = coverage_density, lower.tail = FALSE)
  p_adjusted <- p.adjust(p_value, method = "BH")
  
  theBins$E <- E
  theBins$pvalue <- p_value
  theBins$p_adjusted <- p_adjusted
  theBins$logP <- -log10(p_adjusted)
  theBins$logP[which(is.infinite(theBins$logP))] <- 50
  theBins$logP[which(theBins$logP > 50)] <- 50
  
  return(theBins)
}

plotOneChr <- function(RMRes, genomeSize, chrID, binSize){
  theTypes <- unique(RMRes$V11)
  nt <- length(theTypes)
  
  # plots <- list()
  for(i in 1 : nt){
    message(chrID, " ---- ", theTypes[i])
    oneTypeRM <- RMRes[which(RMRes$V11 == theTypes[i]), ]
    xmins <- oneTypeRM$V6
    xmaxs <- oneTypeRM$V7
    
    theBins <- getBinPvalues(oneTypeRM, genomeSize, binSize)
    
    p <- ggplot() + geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 0, ymax = 3),  fill = "white") +
      geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 1, ymax = 2),  fill = "white", color = "black") + 
      geom_rect(aes(xmin = xmins, xmax = xmaxs, ymin = 1, ymax = 2),  fill = "red") + 
      geom_rect(aes(xmin = theBins$start, xmax = theBins$end, ymin = 0.5, ymax = 0.9), 
                fill = "#2F4F4F", alpha = theBins$logP / 50) +
      ggtitle(paste(chrID, " -- ", theTypes[i], sep = "")) + 
      theme_void()
    print(p)
  }
}

chrs <- unique(M_RepeatMasker$V5)
pdf("./M_repeatMasker_plot_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
 oneRM <- M_RepeatMasker[which(M_RepeatMasker$V5 == oneChr), ]
 oneSize <- Msize$V2[which(Msize$V1 == oneChr)]
 plotOneChr(oneRM, oneSize, oneChr, binSize = 100000)
}
dev.off()

chrs <- unique(W_RepeatMasker$V5)
pdf("./W_repeatMasker_plot_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
  oneRM <- M_RepeatMasker[which(M_RepeatMasker$V5 == oneChr), ]
  oneSize <- Msize$V2[which(Msize$V1 == oneChr)]
  plotOneChr(oneRM, oneSize, oneChr, binSize = 100000)
}
dev.off()








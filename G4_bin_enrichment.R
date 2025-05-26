## repeat masker processing

setwd("~/mode-species/fly/G4/")

library(data.table)
library(ggplot2)
library(cowplot)

rm(list = ls())


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
  
  repeatRange <- IRanges(start = oneTypeRM$V4, end = oneTypeRM$V5)
  xmins <- oneTypeRM$V4
  xmaxs <- oneTypeRM$V5
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


plotOneChr <- function(G4Res, genomeSize, chrID, binSize){
  message(chrID)
  xmins <- G4Res$V4
  xmaxs <- G4Res$V5
  
  theBins <- getBinPvalues(G4Res, genomeSize, binSize)
  
  p <- ggplot() + geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 0, ymax = 3),  fill = "white") +
    geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 1, ymax = 2),  fill = "#F5F5F5", color = "black") + 
    geom_rect(aes(xmin = xmins, xmax = xmaxs, ymin = 1, ymax = 2),  fill = "orange") + 
    geom_rect(aes(xmin = theBins$start, xmax = theBins$end, ymin = 0.5, ymax = 0.9), 
              fill = "#2F4F4F", alpha = theBins$logP / 50) +
    ggtitle(paste(chrID, " -- G4, Score > 100", sep = "")) + 
    theme_void()
  
  print(p)
}

Msize <- fread("../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)
chrs <- paste("chr", c("2L", "2R", "3L", "3R", "4", "X", "Y"), sep = '')
pdf("./M_G4_pqsfinder_plot_100_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
  M_G4_onechr <- fread(paste("./", "M_", oneChr, ".pqsfinder.gff", sep = ""), skip = 3, stringsAsFactors = F, fill = T, sep = "\t")
  M_G4_onechr <- M_G4_onechr[which(M_G4_onechr$V6 > 100), ]
  
  oneSize <- Msize$V2[which(Msize$V1 == oneChr)]
  plotOneChr(M_G4_onechr, oneSize, oneChr, binSize = 100000)
}
dev.off()

Wsize <- fread("../W.np2.manule.fa.fai", header = F, stringsAsFactors = F)
chrs <- paste("chr", c("2L", "2R", "3L", "3R", "4", "X"), sep = '')
pdf("./W_G4_pqsfinder_plot_100_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
  W_G4_onechr <- fread(paste("./", "W_", oneChr, ".pqsfinder.gff", sep = ""), skip = 3, stringsAsFactors = F, fill = T, sep = "\t")
  W_G4_onechr <- W_G4_onechr[which(W_G4_onechr$V6 > 100), ]
  
  oneSize <- Wsize$V2[which(Wsize$V1 == oneChr)]
  plotOneChr(W_G4_onechr, oneSize, oneChr, binSize = 100000)
}
dev.off()

#########################################


######## i-motif prediction, iM-Seeker
#####################################
#####################################
setwd("~/mode-species/fly/G4/")

library(data.table)
library(ggplot2)
library(cowplot)

rm(list = ls())

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
  repeatRange <- IRanges(start = oneTypeRM$beg, end = oneTypeRM$end)
  
  xmins <- oneTypeRM$beg
  xmaxs <- oneTypeRM$end
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


plotOneChr <- function(G4Res, genomeSize, chrID, binSize){
  message(chrID)
  xmins <- G4Res$beg
  xmaxs <- G4Res$end
  
  theBins <- getBinPvalues(G4Res, genomeSize, binSize)
  
  p <- ggplot() + geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 0, ymax = 3),  fill = "white") +
    geom_rect(aes(xmin = 1, xmax = genomeSize, ymin = 1, ymax = 2),  fill = "#F5F5F5", color = "black") + 
    geom_rect(aes(xmin = xmins, xmax = xmaxs, ymin = 1, ymax = 2),  fill = "green", alpha = G4Res$predict_score) + 
    geom_rect(aes(xmin = theBins$start, xmax = theBins$end, ymin = 0.5, ymax = 0.9), 
              fill = "#2F4F4F", alpha = theBins$logP / 50) +
    ggtitle(paste(chrID, " -- iMotif, iM-Seeker", sep = "")) + 
    theme_void()
  
  print(p)
}


Msize <- fread("../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)
chrs <- paste("chr", c("2L", "2R", "3L", "3R", "4", "X", "Y"), sep = '')
M_iMotif <- fread("./M_iMseeker_imotif_prediction.csv", stringsAsFactors = F, fill = T, sep = ",")
M_iMotif <- M_iMotif[which(M_iMotif$predict_score > 0.4)]
pdf("./M_iMotif_plot_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
  M_iMotif_onechr <- M_iMotif[which(M_iMotif$chr == oneChr), ]
  oneSize <- Msize$V2[which(Msize$V1 == oneChr)]
  plotOneChr(M_iMotif_onechr, oneSize, oneChr, binSize = 100000)
}
dev.off()


Wsize <- fread("../W.np2.manule.fa.fai", header = F, stringsAsFactors = F)
chrs <- paste("chr", c("2L", "2R", "3L", "3R", "4", "X"), sep = '')
W_iMotif <- fread("./W_iMseeker-imotif_prediction.csv", stringsAsFactors = F, fill = T, sep = ",")
W_iMotif <- W_iMotif[which(W_iMotif$predict_score > 0.4)]
pdf("./W_iMotif_plot_with_p.pdf", width = 10, height = 3)
for(oneChr in chrs){
  W_iMotif_onechr <- W_iMotif[which(W_iMotif$chr == oneChr), ]
  oneSize <- Wsize$V2[which(Wsize$V1 == oneChr)]
  plotOneChr(W_iMotif_onechr, oneSize, oneChr, binSize = 100000)
}
dev.off()


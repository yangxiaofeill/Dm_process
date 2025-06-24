## repeat masker processing
setwd("~/mode-species/fly/published_genome/chr_genome/")

library(data.table)
library(ggplot2)
library(cowplot)
library(IRanges)
library(dplyr)
library(purrr)

rm(list = ls())


Msize <- fread("../../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)

Sign <- "M.np2.manule.noY.fa"
OtherSpeciesList <- c("Drosophila_mauritiana", "Drosophila_sechellia", "Drosophila_simulans",
                      "Drosophila_teissieri", "Drosophila_yakuba")

pdf("./M_other_species_SYN_plot_2.pdf", width = 10, height = 3)

for(oneSpecies in OtherSpeciesList){
  message(oneSpecies)
  theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
  theChrs <- unique(theSyri$V1)
  theParSyri <- theSyri[1, ]
  theParSyri <- theParSyri[-1, ]
  
  for(oneChr in theChrs){
    if(oneChr == "-"){
      next()
    }
    
    message(oneChr)
    
    oneChrSyri <- theSyri[which(theSyri$V1 == oneChr), ]
    oneAlignTypes <- unique(oneChrSyri$V10)
    oneChrAlign <- oneChrSyri[which(oneChrSyri$V9 %in% oneAlignTypes), ]
    oneChrSV  <- oneChrAlign[which(oneChrAlign$V11 != "SYN" & oneChrAlign$V11 != "INV"), ]

    chrSize <- Msize$V2[which(Msize$V1 == oneChr)]
    
    p <- ggplot() + geom_rect(aes(xmin = 1, xmax = chrSize, ymin = 0, ymax = 3),  fill = "white") +
      geom_rect(aes(xmin = 1, xmax = chrSize, ymin = 1.5, ymax = 2),  fill = "white", color = "black") +
      geom_rect(aes(xmin = as.numeric(oneChrSV$V2), xmax = as.numeric(oneChrSV$V3),
                    ymin = 1.5, ymax = 2),  fill = "#f3796b") +
      
      ggtitle(paste(oneSpecies, " -- ", oneChr, sep = "")) + 
      theme_void()
    print(p)
  }
  
}
dev.off()

######################################
### centromere region identification
######################################

## repeat masker processing
setwd("~/mode-species/fly/published_genome/chr_genome/")
library(data.table)
library(ggplot2)
library(cowplot)
library(IRanges)
library(dplyr)
library(purrr)

rm(list = ls())


Msize <- fread("../../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)

Sign <- "M.np2.manule.noY.fa"
OtherSpeciesList <- c("Drosophila_mauritiana", "Drosophila_sechellia", "Drosophila_simulans",
                      "Drosophila_teissieri", "Drosophila_yakuba")

cenM <- data.frame(chr = "chr2L", start = 21467129, end = 24305334)
cenM <- rbind(cenM, data.frame(chr = "chr2R", start = 1, end = 15609790))
cenM <- rbind(cenM, data.frame(chr = "chr3L", start = 25000000, end = 33891631))
cenM <- rbind(cenM, data.frame(chr = "chr3R", start = 1, end = 1300000))
cenM <- rbind(cenM, data.frame(chr = "chrX", start = 21517507, end = 29816723))

otherCen <- data.frame(species = "", chr = "", start = 0, end = 0, Mstart = 0, Mend = 0)
otherCen <- otherCen[-1, ]

getOneSpeciesOneChrSyn <- function(theSyri, oneChr){
  theChrs <- unique(theSyri$V1)
  theParSyri <- theSyri[1, ]
  theParSyri <- theParSyri[-1, ]
  
  genomeSize <- fread(paste("./", oneSpecies, ".fa.chr_new.fa.fai", sep = ""), header = F, stringsAsFactors = F)
  
  message(oneChr)
  
  otherSpeciesChrSize <- as.numeric(genomeSize$V2[which(genomeSize$V1 == oneChr)])
  
  message(oneSpecies, "---", oneChr, "----", otherSpeciesChrSize)
  oneChrSyri <- theSyri[which(theSyri$V1 == oneChr), ]
  oneAlignTypes <- unique(oneChrSyri$V10)
  oneChrAlign <- oneChrSyri[which(oneChrSyri$V9 %in% oneAlignTypes), ]
  oneChrSYN <- oneChrAlign[which(oneChrAlign$V11 == "SYN" | oneChrAlign$V11 == "INV"), ]
  
  chrSize <- Msize$V2[which(Msize$V1 == oneChr)]
  oneChrCenStart <- cenM$start[which(cenM$chr == oneChr)]
  oneChrCenEnd   <- cenM$end[which(cenM$chr == oneChr)]
  
  return(oneChrSYN)
}

###########################################################################################################################
oneSpecies <- "Drosophila_mauritiana"
############################################################################################################################
theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
oneChr <- "chr2L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21070370 , end = 24228064, Mstart = 22155614, Mend = 24305334))
oneChr <- "chr2R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 6690235, Mstart = 1, Mend = 16067730))
oneChr <- "chr3L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 23270450, end = 26143691, Mstart = 25407101, Mend = 33891632))
oneChr <- "chr3R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 1594648, Mstart = 1, Mend = 1280898))
oneChr <- "chrX"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21310850, end = 23222241, Mstart = 21810467, Mend = 27642780))
###############################################################################################################################

###############################################################################################################################
oneSpecies <- "Drosophila_sechellia"
###############################################################################################################################
theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
oneChr <- "chr2L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21098840, end = 24956976, Mstart = 21571076, Mend = 24305334))
oneChr <- "chr2R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 4587241, Mstart = 1, Mend = 16230764))
oneChr <- "chr3L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 24044062, end = 28131630, Mstart = 25129687, Mend = 33891632))
oneChr <- "chr3R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 3058027, Mstart = 1, Mend = 1080990))
oneChr <- "chrX"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21663433, end = 22892653, Mstart = 21847002, Mend = 28998739))
################################################################################################################################

################################################################################################################################
oneSpecies <- "Drosophila_simulans"
################################################################################################################################
theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
oneChr <- "chr2L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21041174, end = 23857595, Mstart = 21455322, Mend = 24305334))
oneChr <- "chr2R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 5870645, Mstart = 1, Mend = 16404770))
oneChr <- "chr3L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 23379671, end = 23399903, Mstart = 24448128, Mend = 33891632))
oneChr <- "chr3R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 719380, Mstart = 1, Mend = 1080990))
oneChr <- "chrX"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 21017200, end = 22007221, Mstart = 21847005, Mend = 24523608))
##############################################################################################################################


##############################################################################################################################
oneSpecies <- "Drosophila_teissieri"
##############################################################################################################################
theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
oneChr <- "chr2L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 25492899, end = 28277466, Mstart = 23029655, Mend = 24305334))
oneChr <- "chr2R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 3104125, Mstart = 1, Mend = 12597529))
oneChr <- "chr3L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 24401095, end = 24551875, Mstart = 25381600, Mend = 33891632))
oneChr <- "chr3R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 3467571, Mstart = 1, Mend = 1304629))
oneChr <- "chrX"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 20987825, end = 23378501, Mstart = 21851658, Mend = 27529931))
##################################################################################################################################

#################################################################################################################################
oneSpecies <- "Drosophila_yakuba"
##################################################################################################################################
theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
oneChr <- "chr2L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 24253376, end = 29937316, Mstart = 23025747, Mend = 24305334))
oneChr <- "chr2R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 6819501, Mstart = 1, Mend = 18472962))
oneChr <- "chr3L"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 24075227, end = 30213875, Mstart = 25381518, Mend = 33891632))
oneChr <- "chr3R"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 1 , end = 2947335, Mstart = 1, Mend = 1280927))
oneChr <- "chrX"
speChrSYN <- getOneSpeciesOneChrSyn(theSyri, oneChr)
print(speChrSYN)
otherCen <- rbind(otherCen, 
                  data.frame(species = oneSpecies, chr = oneChr, 
                             start = 20780857, end = 23044442, Mstart = 21231986, Mend = 28369404))

######################################################################################################################
write.table(otherCen, file = "./0_Cen_Syn_M_and_other.txt", quote = F, sep = "\t", row.names = F, col.names = F)


###############################################################################################
## programming, does not work
##############################################################################################
for(oneSpecies in OtherSpeciesList){
  message(oneSpecies)
  theSyri <- read.table(paste(Sign, "-", oneSpecies, ".fa.chr_new.fa.syri.out", sep = ''), header = F, stringsAsFactors = F)
  theChrs <- unique(theSyri$V1)
  theParSyri <- theSyri[1, ]
  theParSyri <- theParSyri[-1, ]
  
  genomeSize <- fread(paste("./", oneSpecies, ".fa.chr_new.fa.fai", sep = ""), header = F, stringsAsFactors = F)
  
  for(oneChr in theChrs){
    if(oneChr == "-" | oneChr == "chr4"){
      next()
    }
    message(oneChr)
    
    otherSpeciesChrSize <- as.numeric(genomeSize$V2[which(genomeSize$V1 == oneChr)])
    
    oneChrSyri <- theSyri[which(theSyri$V1 == oneChr), ]
    oneAlignTypes <- unique(oneChrSyri$V10)
    oneChrAlign <- oneChrSyri[which(oneChrSyri$V9 %in% oneAlignTypes), ]
    oneChrSYN <- oneChrAlign[which(oneChrAlign$V11 == "SYN" | oneChrAlign$V11 == "INV"), ]
    
    chrSize <- Msize$V2[which(Msize$V1 == oneChr)]
    oneChrCenStart <- cenM$start[which(cenM$chr == oneChr)]
    oneChrCenEnd   <- cenM$end[which(cenM$chr == oneChr)]
    
    if(oneChr == "chr2L" | oneChr == "chr3L"){
      ## the last syn region with start < M_cen_start
      otherChrCenEnd <- otherSpeciesChrSize
      smallerIndex <- which(as.numeric(oneChrSYN$V2) < oneChrCenStart)
      otherChrCenStart <- as.numeric(oneChrSYN$V7[smallerIndex[length(smallerIndex)]])
      
      MCenSynStart <- as.numeric(oneChrSYN$V2[smallerIndex[length(smallerIndex)]])
      MCenSynEnd   <- as.numeric(Msize$V2[which(Msize$V1 == oneChr)])
    }
    if(oneChr == "chr2R" | oneChr == "chr3R"){
      ## the first syn region with start > M_cen_end
      largerIndex <- which(as.numeric(oneChrSYN$V2) > oneChrCenEnd)
      otherChrCenEnd <- as.numeric(oneChrSYN$V7[largerIndex[1]])
      otherChrCenStart <- 1
      
      MCenSynStart <- 1
      MCenSynEnd   <- as.numeric(oneChrSYN$V2[largerIndex[1]])
    }
    
    if(oneChr == "chrX"){
      smallerIndex <- which(oneChrSYN$V2 < oneChrCenStart)
      otherChrCenStart <- as.numeric(oneChrSYN$V7[smallerIndex[length(smallerIndex)]])
      MCenSynStart <- as.numeric(oneChrSYN$V2[smallerIndex[length(smallerIndex)]])
      
      largerIndex <- which(oneChrSYN$V2 > oneChrCenEnd)
      otherChrCenEnd <- as.numeric(oneChrSYN$V7[largerIndex[1]])
      MCenSynEnd   <- as.numeric(oneChrSYN$V2[largerIndex[1]])
    }
    
    otherCen <- rbind(otherCen, 
                      data.frame(species = oneSpecies, chr = oneChr, 
                                 start = otherChrCenStart, end = otherChrCenEnd, 
                                 Mstart = MCenSynStart, Mend = MCenSynEnd))
  }
}
##############################################################################################

write.table(otherCen, file = "./0_Cen_Syn_M_and_other.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#####################################################################################################################

######################################################################################################################
## the repetative elements for centromere regions
######################################################################################################################
rm(list=ls())
library(data.table)
library(ggplot2)
library(cowplot)
library(IRanges)
library(dplyr)
library(ggpubr)
library(purrr)

setwd("~/mode-species/fly/published_genome/chr_genome/")
Msize <- fread("../../M.np2.manule.fa.fai", header = F, stringsAsFactors = F)
M_RepeatMasker <- fread("../../RepeatMasker/M.np2.manule.fa.out", skip = 3, stringsAsFactors = F, fill = T)

cenPos <- fread("./0_Cen_Syn_M_and_other.txt", header = F, stringsAsFactors = F)
colnames(cenPos) <- c("species", "chr", "species_start", "species_end", "M_start", "M_end")

getFreq <- function(theRepeatMasker, GenomeSize){
  chrs <- unique(theRepeatMasker$V5)
  freqRes <- data.frame(chr = "", rType = "", size = 0, percentage = 0.0)
  freqRes <- freqRes[-1, ]
  for(oneChr in chrs){
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

OtherSpeciesList <- c("Drosophila_mauritiana", "Drosophila_sechellia", "Drosophila_simulans",
                      "Drosophila_teissieri", "Drosophila_yakuba")
#############################################################################################################################################
## get overall frequency of diff rep
#############################################################################################################################################
# for(oneSpe in OtherSpeciesList){
#   message("***************", oneSpe, ":repeatFreq***************")
#   GenomeSize         <- fread(paste(oneSpe, ".fa.chr_new.fa.fai", sep = ""), header = F, stringsAsFactors = F)
#   GenomeRepeatMasker <- fread(paste(oneSpe, ".fa.chr_new.fa.out", sep = ""), skip = 3, stringsAsFactors = F, fill = T)
#   repeatFreq <- getFreq(theRepeatMasker = GenomeRepeatMasker, GenomeSize = GenomeSize)
#   write.table(repeatFreq, file = paste(oneSpe, ".fa.chr_new.fa.out.freq", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)
# }
#############################################################################################################################################

getCenRepFreq<- function(oneTypeRM, cenStart, cenEnd){
  cenRange <- IRanges(start = cenStart, end = cenEnd)
  repeatRange <- IRanges(start = oneTypeRM$V6, end = oneTypeRM$V7)

  # calculate the overlap length
  overlaps <- findOverlaps(cenRange, repeatRange)
  intersect_lengths <- width(pintersect(
    cenRange[queryHits(overlaps)],
    repeatRange[subjectHits(overlaps)]
  ))
  coverage_lengths <- tapply(intersect_lengths, queryHits(overlaps), sum)
  if(length(coverage_lengths) == 0){
    coverage_lengths = 0
  }
  
  return(coverage_lengths)
}

chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
repTypes <- c("LINE/I-Jockey", "LINE/R1", "LTR/Copia", "LTR/Gypsy", "LTR/Pao", "RC/Helitron", "DNA/P")
cenRep <- data.frame(species = "", chr = "", repType = "", cenStart = 0, cenEnd = 0, repContent = 0, repPer = 0, 
                     McenStart = 0, McenEnd = 0, MrepContent = 0, MrepPer = 0)
cenRep <- cenRep[-1, ]

for(oneSpe in OtherSpeciesList){
  message("***************", oneSpe, "***************")
  GenomeSize         <- fread(paste(oneSpe, ".fa.chr_new.fa.fai", sep = ""), header = F, stringsAsFactors = F)
  GenomeRepeatMasker <- fread(paste(oneSpe, ".fa.chr_new.fa.out", sep = ""), skip = 3, stringsAsFactors = F, fill = T)
  for(oneRepType in repTypes){
    for(oneChr in chrs){
      message("***************", oneRepType, "--", oneChr)
      cenIndex <- which(cenPos$species == oneSpe & cenPos$chr == oneChr)
      oneTypeChrRM  <- GenomeRepeatMasker[which(GenomeRepeatMasker$V5 == oneChr & GenomeRepeatMasker$V11 == oneRepType), ]
      oneTypeCenCov <- getCenRepFreq(oneTypeChrRM, cenPos$species_start[cenIndex], cenPos$species_end[cenIndex])
      oneTypeCenPer <- oneTypeCenCov / (cenPos$species_end[cenIndex] - cenPos$species_start[cenIndex] + 1)
      
      MTypeChrRM   <- M_RepeatMasker[which(M_RepeatMasker$V5 == oneChr & M_RepeatMasker$V11 == oneRepType), ]
      MTypeCenCov  <- getCenRepFreq(MTypeChrRM, cenPos$M_start[cenIndex], cenPos$M_end[cenIndex])
      MTypeCenPer  <- MTypeCenCov / (cenPos$M_end[cenIndex] - cenPos$M_start[cenIndex] + 1)
      
      cenRep <- rbind(cenRep,
                      data.frame(species = oneSpe, chr = oneChr, repType = oneRepType, 
                                 cenStart = cenPos$species_start[cenIndex], cenEnd = cenPos$species_end[cenIndex], repContent = oneTypeCenCov, repPer = oneTypeCenPer, 
                                 McenStart = cenPos$M_start[cenIndex], McenEnd = cenPos$M_end[cenIndex], MrepContent = MTypeCenCov, MrepPer = MTypeCenPer))
    }
  }
}

#write.table(cenRep, file = "./0_Cen_Syn_M_and_other_repeat_content.txt", quote = F, sep = "\t", row.names = F, col.names = T)
cenRep <- cenRep[which(cenRep$cenEnd - cenRep$cenStart > 100000), ]
cenRep <- cbind(cenRep, fold = cenRep$repPer / cenRep$MrepPer)

speciesOrder <- c("Drosophila_teissieri", "Drosophila_yakuba", "Drosophila_mauritiana", "Drosophila_sechellia", "Drosophila_simulans")
pdf("./0_Cen_Syn_M_and_other_repeat_content_fold-2.pdf", width = 4, height = 4)
for(oneT in repTypes){
  message("***************", oneT, "***************")
  oneCenRep <- cenRep[which(cenRep$repType == oneT), ]
  
  oneCenRep$species <- factor(oneCenRep$species, levels = speciesOrder)
  
  median_data <- oneCenRep %>% 
    group_by(species) %>%
    summarise(median_value = median(fold))
  
  p <- ggboxplot(oneCenRep, x = "species", y = "fold", color = "species", palette = "nejm", add = "jitter", 
            title = oneT, xlab = "", ylab = "Fold") + rotate_x_text(45) + 
     scale_y_continuous(limits = c(0, max(oneCenRep$fold)), labels = scales::number_format(accuracy = 0.1)) +
    theme(legend.position = 'none') +
    # geom_smooth(data = median_data, 
    #             aes(x = as.numeric(species), y = median_value), 
    #             method = "loess", color = "black", se = F) +
    geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "#53565C", linewidth = 1.0)
  print(p)
}
dev.off()



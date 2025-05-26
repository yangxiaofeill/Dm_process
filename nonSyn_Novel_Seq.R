setwd("~/mode-species/fly/syri/")

library(data.table)

rm(list = ls())

Sign <- "W"

Gsize <- fread(paste("../", Sign, ".np2.manule.fa.fai", sep = ''), header = F, stringsAsFactors = F)
theSyri <- read.table(paste("./", Sign, ".GCF_000001215.4_ref.nucmer.syri.syri.out", sep = ''), header = F, stringsAsFactors = F)
theParSyri <- theSyri[1, ]
theParSyri <- theParSyri[-1, ]
theChrs <- unique(theSyri$V1)
for(oneChr in theChrs){
  message(oneChr)
  oneChrSyri <- theSyri[which(theSyri$V1 == oneChr), ]
  oneParTypes <- unique(oneChrSyri$V10)
  oneChrParSyri <- oneChrSyri[which(oneChrSyri$V9 %in% oneParTypes), ]
  theParSyri <- rbind(theParSyri, oneChrParSyri)
}
write.table(theParSyri, file = paste("./", Sign, ".processed.Syri", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")

# chr4
novelBed <- data.frame(chr = "chr4", start = 1298248, end = 1500479)
# chrX
novelBed <- rbind(novelBed, data.frame(chr = "chrX", start = 25067547, end = 39099338))
# chr3L
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 24356755, end = 27713189))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 30437656, end = 30926382))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 31839081, end = 35593110))
# chr2R 
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 537967, end = 1122714))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 1155615, end = 1892731))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 1926876, end = 2696787))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 3040368, end = 3825498))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 3885746, end = 3962072))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 3985810, end = 4032335))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 4057932, end = 4514373))

novelBed <- novelBed[which(novelBed$end - novelBed$start > 100000), ]
write.table(novelBed, file = paste("./", Sign, ".Novel.Region.bed", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")


Sign <- "M"

Gsize <- fread(paste("../", Sign, ".np2.manule.fa.fai", sep = ''), header = F, stringsAsFactors = F)
theSyri <- read.table(paste("./", Sign, ".GCF_000001215.4_ref.nucmer.syri.syri.out", sep = ''), header = F, stringsAsFactors = F)
theParSyri <- theSyri[1, ]
theParSyri <- theParSyri[-1, ]
theChrs <- unique(theSyri$V1)
for(oneChr in theChrs){
  message(oneChr)
  oneChrSyri <- theSyri[which(theSyri$V1 == oneChr), ]
  oneParTypes <- unique(oneChrSyri$V10)
  oneChrParSyri <- oneChrSyri[which(oneChrSyri$V9 %in% oneParTypes), ]
  theParSyri <- rbind(theParSyri, oneChrParSyri)
}
write.table(theParSyri, file = paste("./", Sign, ".processed.Syri", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")

# chr4
novelBed <- data.frame(chr = "chr4", start = 1307742, end = 1381625)
# chrX
novelBed <- rbind(novelBed, data.frame(chr = "chrX", start = 25345681, end = 35972434))
# chr3R
novelBed <- rbind(novelBed, data.frame(chr = "chr3R", start = 890663, end = 1253438))
novelBed <- rbind(novelBed, data.frame(chr = "chr3R", start = 9603489, end = 9712518))
novelBed <- rbind(novelBed, data.frame(chr = "chr3R", start = 29451835, end = 29504324))
# chr3L
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 7672898, end = 7727152))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 27001316, end = 27060513))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 27312436, end = 27801162))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 28713892, end = 32576780))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 32599760, end = 32781194))
novelBed <- rbind(novelBed, data.frame(chr = "chr3L", start = 33863404, end = 33891632))
# chr2R 
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 5555090, end = 11475447))
novelBed <- rbind(novelBed, data.frame(chr = "chr2R", start = 33179963, end = 33224741))

novelBed <- novelBed[which(novelBed$end - novelBed$start > 100000), ]
write.table(novelBed, file = paste("./", Sign, ".Novel.Region.bed", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")



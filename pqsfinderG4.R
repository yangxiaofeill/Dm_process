library("pqsfinder")
library(rtracklayer)

##SIGN="W"
SIGN <- c("M", "W")
for(oneSign in SIGN){
  message(oneSign)
  fastaFile <- paste("../", oneSign, ".np2.manule.fa", sep = '')
  dnaSeq <- readDNAStringSet(fastaFile)
  chrs <- names(dnaSeq)
  for(i in 1 : length(chrs)){
    message(chrs[i])
    oneseq <- DNAString(as.character(dnaSeq[i]))
    onepqs <- pqsfinder(oneseq, min_score = 25)
    export(onepqs, paste(oneSign, "_", chrs[i], ".pqsfinder.gff", sep = ""), version = "3")
    writeXStringSet(as(onepqs, "DNAStringSet"), file = paste(oneSign, "_", chrs[i], ".pqsfinder.fa", sep = ""), format = "fasta")
  }
}


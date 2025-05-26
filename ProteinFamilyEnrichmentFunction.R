library(data.table)
library(clusterProfiler)
library(ggpubr)

setwd("~/mode-species/fly/Figures/")
### Dm-XJTU-M
rm(list=ls())
interProscanRes <- fread("../geneAnnotation/M_diffType/M_maker.all.protein.iprscan_5.67.final", header = F, stringsAsFactors = F)

novelGff <- fread("../syri/M.Novel.Region.gene.gff", header = F, stringsAsFactors = F)
novelGff <- novelGff[which(novelGff$V3 == "mRNA"), ]
geneIDs <- sapply(strsplit(novelGff$V9, ";"), function(x) x[1])
geneIDs <- sapply(geneIDs, function(x) substr(x, 4, nchar(x)))

interNovelRes <- interProscanRes[which(interProscanRes$V1 %in% geneIDs), ]

freqTerm <- as.data.frame(table(interNovelRes$V6))
freqTerm <- freqTerm[order(freqTerm$Freq, decreasing = T), ]

freqTermTop5 <- freqTerm[1:5, ]
freqTermTop5$Var1 <- factor(freqTermTop5$Var1, levels = freqTermTop5$Var1)

ggbarplot(freqTermTop5, "Var1", "Freq",
          orientation = "horiz",
          fill = "#F39B7E", color = "black",
          label = TRUE, lab.pos = "out", lab.col = "black", 
          ylab = "gene count", xlab = "")


gene2Term <- interProscanRes[, c("V13", "V1")]

enrich_result <- enricher(gene = geneIDs, 
                          TERM2GENE = gene2Term,
                          TERM2NAME = NA,
                          pAdjustMethod = "BH")
sigPFam <- enrich_result@result
sigPFam <- sigPFam[which(sigPFam$p.adjust < 0.05), ]
log10Padj <- -log10(sigPFam$p.adjust)
sigPFam$log10Padj <- log10Padj

ggbarplot(sigPFam, "ID", "log10Padj",
          orientation = "horiz",
          fill = "#F39B7E", color = "black",
          label = FALSE, lab.pos = "out", lab.col = "black", 
          ylab = "-log10(p.adj)", xlab = "") +
  geom_hline(
    yintercept = -log10(0.05), 
    color = "black",
    linetype = "dashed",
    linewidth = 1
  )




### Dm-XJTU-F
rm(list=ls())
interProscanRes <- fread("../geneAnnotation/W_diffType/W_maker.all.protein.iprscan_5.67.final", header = F, stringsAsFactors = F)

novelGff <- fread("../syri/W.Novel.Region.gene.gff", header = F, stringsAsFactors = F)
novelGff <- novelGff[which(novelGff$V3 == "mRNA"), ]
geneIDs <- sapply(strsplit(novelGff$V9, ";"), function(x) x[1])
geneIDs <- sapply(geneIDs, function(x) substr(x, 4, nchar(x)))

interNovelRes <- interProscanRes[which(interProscanRes$V1 %in% geneIDs), ]
freqTerm <- as.data.frame(table(interNovelRes$V6))
freqTerm <- freqTerm[order(freqTerm$Freq, decreasing = T), ]

freqTermTop5 <- freqTerm[1:5, ]
freqTermTop5$Var1 <- factor(freqTermTop5$Var1, levels = freqTermTop5$Var1)

ggbarplot(freqTermTop5, "Var1", "Freq",
          orientation = "horiz",
          fill = "#009F86", color = "black",
          label = TRUE, lab.pos = "out", lab.col = "black", 
          ylab = "gene count", xlab = "")


gene2Term <- interProscanRes[, c("V6", "V1")]

enrich_result <- enricher(gene = geneIDs, 
                          TERM2GENE = gene2Term,
                          TERM2NAME = NA,
                          pAdjustMethod = "BH")
sigPFam <- enrich_result@result
sigPFam <- sigPFam[which(sigPFam$p.adjust < 0.05), ]
log10Padj <- -log10(sigPFam$p.adjust)
sigPFam$log10Padj <- log10Padj

ggbarplot(sigPFam, "ID", "log10Padj",
          orientation = "horiz",
          fill = "#009F86", color = "black",
          label = FALSE, lab.pos = "out", lab.col = "black", 
          ylab = "-log10(p.adj)", xlab = "") + 
  geom_hline(
    yintercept = -log10(0.05), 
    color = "black",
    linetype = "dashed",
    linewidth = 1
  )



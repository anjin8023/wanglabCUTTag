setwd("/home/liukuai/Documents/project/shenx/Cut_Tag_anti-CFP1")
library(GenomicRanges)
library(dplyr)
library(stringr)
library(ChIPseeker)
library(org.Mm.eg.db)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

files = list.files("different_binding_site/macs2", full.names = T)

for (i in files){
  filename = basename(i)
  peakRes = readPeakFile(i)
  keepchr = !grepl("[.]", seqlevels(peakRes))
  seqlevels(peakRes, pruning.mode="coarse") <- seqlevels(peakRes)[keepchr]
  peakAnno = annotatePeak(peakRes, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno2 = as.data.frame(peakAnno)
  
  pdf(paste0(filename,".pdf"))
  plotAnnoPie(peakAnno)
  dev.off()
  #xlsx::write.xlsx(peakAnno2, file = paste0(filename, "_peaks.xlsx"))
}




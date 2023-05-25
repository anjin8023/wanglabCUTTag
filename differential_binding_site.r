setwd("~/Documents/project/shenx/Cut_Tag_anti-CFP1/")

library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(chromVAR) 
library(DESeq2) 
library(ggpubr) 
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene


mPeak = GRanges()
peaksDir='different_binding_site/macs2/'
histL = list.files(peaksDir)
histL = basename(histL)

## Create a master peak list merging all the peaks called for each sample
for(hist in histL){
  peakRes = read.table(paste0(peaksDir , hist), header = FALSE, fill = TRUE)
  peakRes = peakRes[grep("chr", peakRes[,1]), ]
  mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}

masterPeak = reduce(mPeak)

peakAnno = annotatePeak(masterPeak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno = as.data.frame(peakAnno)
rownames(peakAnno) <- paste(peakAnno$seqnames, peakAnno$start, peakAnno$end, sep = "_")
## Get the fragment counts for each peak in the master peak list

bamDir = "different_binding_site/bam/"
histL2 = list.files(bamDir, pattern = '*.bam$')
histL2 = basename(histL2)

countMat = matrix(NA, length(masterPeak), length(histL2))

## overlap with bam file to get count
i = 1
for(hist in histL2){
  bamFile = paste0(bamDir, "/", hist)
  fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i] = counts(fragment_counts)[,1]
  i = i + 1
}

colnames(countMat) = histL2
colnames(countMat) = gsub('.bam','',colnames(countMat))

index = as.data.frame(masterPeak)
rownames(countMat) = paste(index$seqnames, index$start, index$end, sep = "_")
countMat = countMat[! grepl('^chrMT', rownames(countMat)), ]

countMat_ccr6 <- countMat[, grep("CCR6", colnames(countMat))]
countMat_dn <- countMat[, grep("DN", colnames(countMat))]
countMat_ncr <- countMat[, grep("NCR|NKp46", colnames(countMat))]

##------------------------------##

dat <- countMat_ncr

selectR = which(rowSums(dat) > 5) ## remove low count genes
dataS = dat[selectR,]
condition = factor(rep(c('KO','WT'), times=2), levels = c('WT',"KO"))

dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
countMatDiff = cbind(dataS, normDDS, res)
table(countMatDiff$pvalue < 0.05)
head(countMatDiff)

countMatDiff <- merge(countMatDiff, peakAnno, by=0)


write.csv(peakAnno, file="different_binding_site/results/peakAnno.csv", row.names = F, quote = F)
write.csv(countMatDiff, file="different_binding_site/results/ccr6_countMatDiff.csv", row.names = F, quote = F)
write.csv(countMatDiff, file="different_binding_site/results/dn_countMatDiff.csv", row.names = F, quote = F)
write.csv(countMatDiff, file="different_binding_site/results/ncr_countMatDiff.csv", row.names = F, quote = F)



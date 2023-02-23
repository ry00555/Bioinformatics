setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
```{r calculate TPM}
library("scater")
library("dplyr")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
BiocManager::install("DESeq2")
BiocManager::install("scales")
BiocManager::install("d3heatmap")
library("d3heatmap")
BiocManager::install("IHW")
library("ggplot2")
library("scales")
library("readr")

ReportersSubJGIAllcountsKOoI <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.txt", row.names = 1, check.names = FALSE, header=TRUE, sep="\t")
ReportersSubJGIAllcountsKOoImatrix <- data.matrix(ReportersSubJGIAllcountsKOoI)
ReportersSubJGIAllcountsKOoIcountsOnly <- data.matrix(ReportersSubJGIAllcountsKOoImatrix[ ,6:ncol(ReportersSubJGIAllcountsKOoImatrix)])


#calculateTPM
ReportersSubJGIAllcountsKOoITPMplus <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly, lengths = ReportersSubJGIAllcountsKOoI[,5])+1
##convert to matrix again
ReportersSubJGIAllcountsKOoITPMplus  <- data.matrix(ReportersSubJGIAllcountsKOoITPMplus)


countData <- read.csv('/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv', stringsAsFactors=FALSE, header = TRUE, row.names=1,  sep = ",")
head(countDataMatrixplus)
countDataMatrixplus <- as.matrix(countData)+1
countDatas_numeric <- matrix(as.numeric(countDataMatrixplus), ncol = 85)
dimnames(countDatas_numeric) <- list(rownames(countDataMatrixplus), colnames(countDataMatrixplus))
CountdataITPMplus <- calculateTPM(countDataMatrixplus, lengths = countData[,5])
##convert to matrix again
CountdataITPMplus  <- data.matrix(CountdataITPMplus)

metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
rownames(metadata) <- make.names(metadata[,1], unique = TRUE)


Ordered_KO_data <- cbind(ReportersSubJGIAllcountsKOoITPMplus[,1:3],ReportersSubJGIAllcountsKOoITPMplus[,4:5],ReportersSubJGIAllcountsKOoITPMplus[,6:8],
                         ReportersSubJGIAllcountsKOoITPMplus[,9:11],ReportersSubJGIAllcountsKOoITPMplus[,12:14],ReportersSubJGIAllcountsKOoITPMplus[,15:17],
                         ReportersSubJGIAllcountsKOoITPMplus[,18:20],ReportersSubJGIAllcountsKOoITPMplus[,21:23],ReportersSubJGIAllcountsKOoITPMplus[,24:26],
                         ReportersSubJGIAllcountsKOoITPMplus[,27:29],ReportersSubJGIAllcountsKOoITPMplus[,30:32],ReportersSubJGIAllcountsKOoITPMplus[,33:35],
                         ReportersSubJGIAllcountsKOoITPMplus[,36:38],ReportersSubJGIAllcountsKOoITPMplus[,39:41],ReportersSubJGIAllcountsKOoITPMplus[,42:44],
                         ReportersSubJGIAllcountsKOoITPMplus[,45:47],ReportersSubJGIAllcountsKOoITPMplus[,48:50],ReportersSubJGIAllcountsKOoITPMplus[,51:53],
                         ReportersSubJGIAllcountsKOoITPMplus[,54:56],ReportersSubJGIAllcountsKOoITPMplus[,57:58],ReportersSubJGIAllcountsKOoITPMplus[,59:61],
                         ReportersSubJGIAllcountsKOoITPMplus[,62:64],ReportersSubJGIAllcountsKOoITPMplus[,65:67],ReportersSubJGIAllcountsKOoITPMplus[,68:70],
                         ReportersSubJGIAllcountsKOoITPMplus[,71:73],ReportersSubJGIAllcountsKOoITPMplus[,74:76],ReportersSubJGIAllcountsKOoITPMplus[,77:79],
                         ReportersSubJGIAllcountsKOoITPMplus[,80:82],ReportersSubJGIAllcountsKOoITPMplus[,83:85])


Averaged_Orderd_KO_data <- cbind(rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,1:3], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,4:5], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,6:8], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,9:11], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,12:14], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,15:17], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,18:20], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,21:23], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,24:26], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,27:29], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,30:32], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,33:35], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,36:38], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,39:41], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,42:44], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,45:47], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,48:50], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,51:53], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,54:56], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,57:58], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,59:61], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,62:64], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,65:67], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,68:70], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,71:72], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,74:76], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,77:79], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,80:82], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPMplus[,83:85], na.rm = TRUE))

averageRowIDs =c("set7", "isw", "suz12", "eed", "crf4.3", "ncu00548", "ncu06787", "ncu06788", "H3K4R", "ncu06788", "ncu04017", "ncu00548", "cac3", "cac2", "ncu00423", "ncu00423", "ncu00548", "set2", "ncu09120", "ncu03481", "ncu03481", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017", "rtt109", "suz12", "wildtype")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

dds <- DESeqDataSetFromMatrix(countData = ReportersSubJGIAllcountsKOoIcountsOnly, 
                              colData = metadata, 
                              design= ~ condition)
dds$condition <- relevel(dds$condition, ref = "wildtype")
dds <- DESeq(dds)


normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)
res <- results(dds)
write.csv(as.data.frame(res[order(res$padj),] ), file="res.csv")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

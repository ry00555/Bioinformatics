setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI")

#Calculate TPM
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
```{r calculate TPM}
library("scater")
library("dplyr")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("ggplot2")
library("scales")
library("readr")


ReportersSubJGIAllcountsKOoI <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersSubJGIAllcountsKOoImatrix <- data.matrix(ReportersSubJGIAllcountsKOoI)
ReportersSubJGIAllcountsKOoIcountsOnly <- data.matrix(ReportersSubJGIAllcountsKOoImatrix[ ,6:ncol(ReportersSubJGIAllcountsKOoImatrix)])

#calculateTPM
ReportersSubJGIAllcountsKOoITPMplus <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly, lengths = ReportersSubJGIAllcountsKOoI[,5])
##convert to matrix again
ReportersSubJGIAllcountsKOoITPMplus  <- data.matrix(ReportersSubJGIAllcountsKOoITPMplus)

ReportersnonClusteredKOoI<- pheatmap(ReportersSubJGIAllcountsKOoITPMplus, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(200), cellwidth = 20, cellheight = 10, legend=T, show_rownames=T, show_colnames=T, fontsize_col=8, fontsize_row = 7, treeheight_row=0, treeheight_col=20, height = 1.5, width = 1.5)

samplesname <- colnames(ReportersSubJGIAllcountsKOoITPMplus)
write.table(samplesname, file="samplesname.txt", sep="/t")

## Importing featurecounts matrix
Read in the featurecounts output, convert to a matrix, and remove unneccesary columns to include ONLY the count numbers, then convert columns to numeric.
```{r cts}

countData <- read.csv('/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv', stringsAsFactors=FALSE, header = TRUE, row.names=1,  sep = ",")
head(countDataMatrixplus)
countDataMatrixplus <- as.matrix(countData)+1
countDatas_numeric <- matrix(as.numeric(countDataMatrixplus), ncol = 85)
dimnames(countDatas_numeric) <- list(rownames(countDataMatrixplus), colnames(countDataMatrixplus))
CountdataITPMplus <- calculateTPM(countDataMatrixplus, lengths = countData[,5])
##convert to matrix again
CountdataITPMplus  <- data.matrix(CountdataITPMplus)



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
TPMmean <- colMeans(Averaged_Orderd_KO_data)
TPMMedian <- colMedians(Averaged_Orderd_KO_data)
TPMstdev <- apply(Averaged_Orderd_KO_data,2,sd)
w <- cbind(TPMmean,TPMMedian, TPMstdev)

ReportersClusteredKOoI<- pheatmap(Averaged_Orderd_KO_data, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = 10, scale="row", cellheight = 10,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, fontsize_row = 7, treeheight_row=0, treeheight_col=20, height = 8.5, width = 7)


colnames(Ordered_KO_data) <- rownames(metadata)
all((rownames(metadata)) == colnames(Ordered_KO_data))
dds2 <- DESeqDataSetFromMatrix(countData = round(Ordered_KO_data), colData = metadata, design = ~condition)
dds2$condition <- relevel(dds2$condition, ref = 'wildtype')
dds3 <- DESeq(dds2)
res3 <- results(dds3)
resultsNames(dds3)
plotMA(dds3, ylim=c(-2,2))

````````````````````````DESeq Analysis````````````````````````````````````````````````


metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
rownames(metadata) <- make.names(metadata[,1], unique = TRUE)
str(metadata)

all(rownames(metadata) == colnames(CountdataITPMplus))

dds <- DESeqDataSetFromMatrix(countData = round(CountdataITPMplus), 
                              colData = metadata, 
                              design= ~ condition)

dds$condition <- relevel(dds$condition, ref = "wildtype")
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
write.csv(as.data.frame(res), file="res.csv")

plotCounts(Averaged_Orderd_KO_data, gene="NCU07149", intgroup = colnames(Averaged_Orderd_KO_data))
ggplot(d, aes(x=condition, y = count)) + geom_point(position=position_jitter(w=0.1, h=0))
plotMA(res, ylim=c(-2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
write.csv(as.data.frame(res[order(res$log2FoldChange),] ), file="res.csv")

d <- plotCounts(dds3, gene="NCU07149", intgroup="condition", returnData = TRUE)
plotCounts(dds, gene="NCU06889", intgroup="condition")
plotCounts(dds3, gene="NCU08085", intgroup="condition")
plotCounts(dds, gene="NCU09953", intgroup="condition")
plotCounts(dds, gene="NCU05106", intgroup="condition")
plotCounts(dds, gene="NCU05107", intgroup="condition")


,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,Feb 22 2023,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
ReportersSubJGIAllcountsKOoI <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.txt", row.names = 1, check.names = FALSE, header=TRUE, sep="\t")
ReportersSubJGIAllcountsKOoImatrix <- data.matrix(ReportersSubJGIAllcountsKOoI)
ReportersSubJGIAllcountsKOoIcountsOnly <- data.matrix(ReportersSubJGIAllcountsKOoImatrix[ ,6:ncol(ReportersSubJGIAllcountsKOoImatrix)])
ReporterGenes <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/ReporterGenes.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")

#calculateTPM
ReportersSubJGIAllcountsKOoITPMplus <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly, lengths = ReportersSubJGIAllcountsKOoI[,5])+1
##convert to matrix again
ReportersSubJGIAllcountsKOoITPMplus  <- data.matrix(ReportersSubJGIAllcountsKOoITPMplus)

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


breaks <- c(seq(0.1, 10, by=0.1))
heatmap4<- pheatmap(Averaged_Orderd_KO_data, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100), breaks = breaks, scale ="row", cellwidth = 10, cellheight = 10,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=20, height = 5, width = 5)
heatmap_plot4<-heatmap4[[4]]
barp <- barplot(Averaged_Orderd_KO_data[1,],las=3, main="NCU05106",col = brewer.pal(23, name="Set3"), ylab= "Average TPM")
barp <- barplot(Averaged_Orderd_KO_data[2,],las=3, main="NCU05107",col = brewer.pal(23, name="Set3"), ylab= "Average TPM")
barp <- barplot(Averaged_Orderd_KO_data[4,],las=3, main="NCU07149",col = brewer.pal(23, name="Set3"), ylab= "Average TPM")



TPMmean <- colMeans(Averaged_Orderd_KO_data)
TPMMedian <- colMedians(Averaged_Orderd_KO_data)
TPMstdev <- apply(Averaged_Orderd_KO_data,2,sd)
w <- cbind(TPMmean,TPMMedian, TPMstdev)

metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
rownames(metadata) <- make.names(metadata[,1], unique = TRUE)
all(rownames(metadata) <- colnames(ReportersSubJGIAllcountsKOoIcountsOnly))
dds <- DESeqDataSetFromMatrix(countData = ReportersSubJGIAllcountsKOoIcountsOnly, 
                              colData = metadata, 
                              design= ~ condition)
````````````````````````````````````````````````````````````````````````````````````````````````````````
ReportersSubJGIAllcountsKOoI <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersSubJGIAllcountsKOoImatrix <- data.matrix(ReportersSubJGIAllcountsKOoI)
ReportersSubJGIAllcountsKOoIcountsOnly <- data.matrix(ReportersSubJGIAllcountsKOoImatrix[ ,6:ncol(ReportersSubJGIAllcountsKOoImatrix)])

#calculateTPM
ReportersSubJGIAllcountsKOoITPMplus <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly, lengths = ReportersSubJGIAllcountsKOoI[,5])
##convert to matrix again
ReportersSubJGIAllcountsKOoITPMplus  <- data.matrix(ReportersSubJGIAllcountsKOoITPMplus)

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

metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
rownames(metadata) <- make.names(metadata[,1], unique = TRUE)
colnames(Ordered_KO_data) <- rownames(metadata)
all((rownames(metadata)) == colnames(Ordered_KO_data))

dds2 <- DESeqDataSetFromMatrix(countData = round(Ordered_KO_data), colData = metadata, design = ~condition)
dds2$condition <- relevel(dds2$condition, ref = "wildtype")
dds3 <- DESeq(dds2)
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/dds3counts")

res_isw <- results(dds3, alpha=alpha, contrast=c("condition", "wildtype", "isw"))
res_set7 <- results(dds3, alpha=alpha, contrast=c("condition", "set7", "wildtype"))
res_rtt109 <- results(dds3, alpha=alpha, contrast=c("condition", "rtt109", "wildtype"))
res_cac3 <- results(dds3, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_H3K4R <- results(dds3, alpha=alpha, contrast=c("condition", "H3K4R", "wildtype"))
res_ncu06788 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu06788", "wildtype"))
res_ncu04017 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu04017", "wildtype"))
res_ncu00548 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu00548", "wildtype"))
res_suz12 <- results(dds3, alpha=alpha, contrast=c("condition", "suz12", "wildtype"))
res_eed <- results(dds3, alpha=alpha, contrast=c("condition", "eed", "wildtype"))
res_crf4.3 <- results(dds3, alpha=alpha, contrast=c("condition", "cac2", "wildtype"))
res_cac2 <- results(dds3, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_ncu00423 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu00423", "wildtype"))
res_ncu06787 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu06787", "wildtype"))
res_set2 <- results(dds3, alpha=alpha, contrast=c("condition", "set2", "wildtype"))
res_ncu09120 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu09120", "wildtype"))
res_ncu03481 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu03481", "wildtype"))
res_ncu02695 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu02695", "wildtype"))
res_ncu03461 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu03461", "wildtype"))
res_nst1 <- results(dds3, alpha=alpha, contrast=c("condition", "nst1", "wildtype"))

write.csv(as.data.frame(res_isw), file = "res_isw.csv")
write.csv(as.data.frame(res_set7), file = "res_set7.csv")
write.csv(as.data.frame(res_rtt109), file = "res_rtt109.csv")
write.csv(as.data.frame(res_cac3), file = "res_cac3.csv")
write.csv(as.data.frame(res_H3K4R), file = "res_H3K4R.csv")
write.csv(as.data.frame(res_ncu06788), file = "res_ncu06788.csv")
write.csv(as.data.frame(res_ncu04017), file = "res_ncu04017.csv")
write.csv(as.data.frame(res_ncu00548), file = "res_ncu00548.csv")
write.csv(as.data.frame(res_suz12), file = "res_suz12.csv")
write.csv(as.data.frame(res_eed), file = "res_eed.csv")
write.csv(as.data.frame(res_crf4.3), file = "res_crf4.3.csv")
write.csv(as.data.frame(res_cac2), file = "res_cac2.csv")
write.csv(as.data.frame(res_ncu00423), file = "res_ncu00423.csv")
write.csv(as.data.frame(res_ncu06787), file = "res_ncu06787.csv")
write.csv(as.data.frame(res_set2), file = "res_set2.csv")
write.csv(as.data.frame(res_ncu09120), file = "res_ncu09120.csv")
write.csv(as.data.frame(res_ncu03481), file = "res_ncu03481.csv")
write.csv(as.data.frame(res_ncu02695), file = "res_ncu02695.csv")
write.csv(as.data.frame(res_ncu03461), file = "res_ncu03461.csv")
write.csv(as.data.frame(res_nst1), file = "res_nst1.csv")

list_of_sigfiles <- list.files(path = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/dds3counts",
                               recursive = TRUE,
                               pattern = ".csv$")
Reporters <- readr::read_csv(list_of_sigfiles, id = "file_name")
l2fc_Reporters <- data.frame(pivot_wider(data = Reporters, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))

rownames(l2fc_Reporters) <- l2fc_Reporters[,1]
l2fc_Reporters <- l2fc_Reporters[,-1]
write.csv(as.data.frame(l2fc_Reporters), file="l2fc_Reporters.csv")
ReportersL2FC <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/dds3counts/l2fc_Reportersdds3.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersL2FCmatrix <- data.matrix(ReportersL2FC)
heatmap4<- pheatmap(ReportersL2FCmatrix, breaks=breaks, color = colorRampPalette((brewer.pal(n = 3, name="OrRd")))(256), scale ="none", cellwidth = 11, cellheight = 15,  cluster_rows =F, cluster_cols = T, clustering_method="complete", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=12, fontsize_row=12, treeheight_row=0, treeheight_col=20, height = 10, width = 10) 

plotCounts(l2fc_Reporters, gene="NCU07149", intgroup="condition", las=1)



summary(dds3)
resultsNames(dds3)
breaks = c
heatmap4<- pheatmap(Averaged_Orderd_KO_data,color = colorRampPalette((brewer.pal(n = 9, name="OrRd")))(256), scale ="none", cellwidth = 11, cellheight = 15,  cluster_rows =F, cluster_cols = T, clustering_method="complete", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=12, fontsize_row=12, treeheight_row=0, treeheight_col=20, height = 10, width = 10) 

write.table(Averaged_Orderd_KO_data, file = "Average Reporter Gene TPM KOoI.txt", col.names = TRUE,
            row.names = TRUE, sep = "\t")
barp <- barplot(Averaged_Orderd_KO_data[6,], main="NCU09953", las = 3,
       ylab="Average TPM", col = brewer.pal(23, name="Set3"))+ text(barp, Averaged_Orderd_KO_data + 0.5, labels=Averaged_Orderd_KO_data )
boxp <- boxplot(Averaged_Orderd_KO_data, las=3, xlab="Strain", ylab="AverageTPM")
head(Averaged_Orderd_KO_data)

````````````````````````DESeq Analysis V2 ````````````````````````````````````````````````
ReportersSubJGIAllcountsKOoI2 <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersSubJGIAllcountsKOoImatrix2 <- data.matrix(ReportersSubJGIAllcountsKOoI2)
ReportersSubJGIAllcountsKOoIcountsOnly2plus <- data.matrix(ReportersSubJGIAllcountsKOoImatrix2[ ,6:ncol(ReportersSubJGIAllcountsKOoI2)])+1

#calculateTPM
ReportersSubJGIAllcountsKOoITPM <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly2, lengths = ReportersSubJGIAllcountsKOoI2[,5])
##convert to matrix again
ReportersSubJGIAllcountsKOoITPM  <- data.matrix(ReportersSubJGIAllcountsKOoITPM)

Averaged_Orderd_KO_TPM <- cbind(rowMeans(ReportersSubJGIAllcountsKOoITPM[,1:3], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,4:5], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,6:8], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,9:11], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,12:14], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,15:17], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,18:20], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,21:23], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,24:26], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,27:29], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,30:32], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,33:35], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,36:38], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,39:41], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,42:44], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,45:47], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,48:50], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,51:53], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,54:56], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,57:58], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,59:61], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,62:64], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,65:67], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,68:70], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,71:72], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,74:76], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,77:79], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,80:82], na.rm = TRUE),
                                 rowMeans(ReportersSubJGIAllcountsKOoITPM[,83:85], na.rm = TRUE))

averageRowIDs =c("set7", "isw", "suz12", "eed", "crf4.3", "ncu00548", "ncu06787", "ncu06788", "H3K4R", "ncu06788", "ncu04017", "ncu00548", "cac3", "cac2", "ncu00423", "ncu00423", "ncu00548", "set2", "ncu09120", "ncu03481", "ncu03481", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017", "rtt109", "suz12", "wildtype")
colnames(Averaged_Orderd_KO_TPM) <- rownames(metadata_averaged)
metadata_averaged <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_reporter_average.csv", check.names = FALSE, header=TRUE, sep=",")
rownames(metadata_averaged) <- make.names(metadata_averaged[,1], unique = TRUE)
all((rownames(metadata_averaged)) == colnames(Averaged_Orderd_KO_TPM))

ddsAVG <- DESeqDataSetFromMatrix(countData = round(Averaged_Orderd_KO_TPM), colData = metadata_averaged, design = ~condition)
ddsAVG$condition <- relevel(ddsAVG$condition, ref = "wildtype")


ddsAVG1 <- DESeq(ddsAVG)
summary(ddsAVG)
resultsNames(dds3)
breaks = c
heatmap4<- pheatmap(Averaged_Orderd_KO_TPM,color = colorRampPalette((brewer.pal(n = 9, name="OrRd")))(256), scale ="none", cellwidth = 11, cellheight = 15,  cluster_rows =F, cluster_cols = T, clustering_method="complete", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=12, fontsize_row=12, treeheight_row=0, treeheight_col=20, height = 10, width = 10) 

colnames(ReportersSubJGIAllcountsKOoIcountsOnly2plus) <- rownames(metadata)
colnames(ReportersSubJGIAllcountsKOoIcountsOnly2plus) == rownames(metadata)
ddsAVG <- DESeqDataSetFromMatrix(countData = round(ReportersSubJGIAllcountsKOoIcountsOnly2plus), colData = metadata, design = ~condition)
ddsAVG$condition <- relevel(ddsAVG$condition, ref = "wildtype")

ddsAVG1 <- DESeq(ddsAVG)
ddsAVG1$condition
ddsAVG1$type <- relevel(ddsAVG1$type, "control")
summary(ddsAVG1)
ddsAVG1name <- write.table(resultsNames(ddsAVG1))
plotCounts(ddsAVG1, gene="NCU07149", intgroup="condition", las=1)

res <- results(ddsAVG1)
head(res)
res
write.csv(as.data.frame(res), file="ddsAVG.csv")
alpha = 0.05
res_isw <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "isw", "wildtype"))
res_set7 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "set7", "wildtype"))
res_rtt109 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "rtt109", "wildtype"))
res_cac3 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_H3K4R <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "H3K4R", "wildtype"))
res_ncu06788 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu06788", "wildtype"))
res_ncu04017 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu04017", "wildtype"))
res_ncu00548 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu00548", "wildtype"))
res_suz12 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "suz12", "wildtype"))
res_eed <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "eed", "wildtype"))
res_crf4.3 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "cac2", "wildtype"))
res_cac2 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_ncu00423 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu00423", "wildtype"))
res_ncu06787 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu06787", "wildtype"))
res_set2 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "set2", "wildtype"))
res_ncu09120 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu09120", "wildtype"))
res_ncu03481 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu03481", "wildtype"))
res_ncu02695 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu02695", "wildtype"))
res_ncu03461 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "ncu03461", "wildtype"))
res_nst1 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "nst1", "wildtype"))

write.csv(as.data.frame(res_isw), file = "res_isw.csv")
write.csv(as.data.frame(res_set7), file = "res_set7.csv")
write.csv(as.data.frame(res_rtt109), file = "res_rtt109.csv")
write.csv(as.data.frame(res_cac3), file = "res_cac3.csv")
write.csv(as.data.frame(res_H3K4R), file = "res_H3K4R.csv")
write.csv(as.data.frame(res_ncu06788), file = "res_ncu06788.csv")
write.csv(as.data.frame(res_ncu04017), file = "res_ncu04017.csv")
write.csv(as.data.frame(res_ncu00548), file = "res_ncu00548.csv")
write.csv(as.data.frame(res_suz12), file = "res_suz12.csv")
write.csv(as.data.frame(res_eed), file = "res_eed.csv")
write.csv(as.data.frame(res_crf4.3), file = "res_crf4.3.csv")
write.csv(as.data.frame(res_cac2), file = "res_cac2.csv")
write.csv(as.data.frame(res_ncu00423), file = "res_ncu00423.csv")
write.csv(as.data.frame(res_ncu06787), file = "res_ncu06787.csv")
write.csv(as.data.frame(res_set2), file = "res_set2.csv")
write.csv(as.data.frame(res_ncu09120), file = "res_ncu09120.csv")
write.csv(as.data.frame(res_ncu03481), file = "res_ncu03481.csv")
write.csv(as.data.frame(res_ncu02695), file = "res_ncu02695.csv")
write.csv(as.data.frame(res_ncu03461), file = "res_ncu03461.csv")
write.csv(as.data.frame(res_nst1), file = "res_nst1.csv")

list_of_sigfiles <- list.files(path = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/ddsAVD1res",
                               recursive = TRUE,
                               pattern = ".csv$")
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/ddsAVD1res")
Reporters <- readr::read_csv(list_of_sigfiles, id = "file_name")
l2fc_Reporters <- data.frame(pivot_wider(data = Reporters, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
plotCounts(l2fc_Reporters, gene="NCU07149", intgroup="condition", las=1)
rownames(l2fc_Reporters) <- l2fc_Reporters[,1]
l2fc_Reporters <- l2fc_Reporters[,-1]
write.csv(as.data.frame(l2fc_Reporters), file="l2fc_Reporters.csv")
ReportersL2FC <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/ddsAVD1res/l2fc_ReportersRENAMED.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersL2FCmatrix <- data.matrix(ReportersL2FC)

breaks <- c(seq(0, 10, by=1))
heatmap4<- pheatmap(ReportersL2FCmatrix, breaks=breaks, color = colorRampPalette((brewer.pal(n = 3, name="OrRd")))(256), scale ="none", cellwidth = 11, cellheight = 15,  cluster_rows =F, cluster_cols = T, clustering_method="complete", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=12, fontsize_row=12, treeheight_row=0, treeheight_col=20, height = 10, width = 10) 
barplot(ReportersL2FCmatrix[3,], las=3)

barp <- barplot(Averaged_Orderd_KO_data[,], main="NCU09953", las = 3,
                ylab="Average TPM", col = brewer.pal(23, name="Set3"))+ text(barp, Averaged_Orderd_KO_data + 0.5, labels=Averaged_Orderd_KO_data )

,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,Feb 22 2023,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
ReportersSubJGIAllcountsKOoI <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersSubJGIAllcountsKOoIMATRIX <- as.matrix(ReportersSubJGIAllcountsKOoI)
ReportersSubJGIAllcountsKOoIcountsOnly <- data.matrix(ReportersSubJGIAllcountsKOoIMATRIX[ ,6:ncol(ReportersSubJGIAllcountsKOoIMATRIX)])
metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
averageRowIDs =c("set7", "isw", "suz12", "eed", "crf4.3", "ncu00548", "ncu06787", "ncu06788", "H3K4R", "ncu06788", "ncu04017", "ncu00548", "cac3", "cac2", "ncu00423", "ncu00423", "ncu00548", "set2", "ncu09120", "ncu03481", "ncu03481", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017", "rtt109", "suz12", "wildtype")
res_isw <- results(ddsAVG1, alpha=alpha, contrast=c("condition", "isw", "wildtype"))

ddsAVG <- DESeqDataSetFromMatrix(countData = ReportersSubJGIAllcountsKOoIcountsOnly, 
                              colData = metadata, 
                              design= ~ condition)
ddsAVG$condition <- relevel(ddsAVG$condition, ref = "wildtype")
dds3 <- DESeq(ddsAVG)


deseq2_results_writer <- function(strain_name, deseq2_obj, output_dir){
  
  mutant_result_1 <- results(ddsAVG1, alpha=alpha, contrast=c("condition", strain_name, "wildtype"))
  
}

list_of_sigfiles <- list.files(path = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/ddsAVD1res",
                               recursive = TRUE,
                               pattern = ".csv$")
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/ddsAVD1res")

results(ddsAVG1, alpha=alpha, contrast=c("condition", "isw", "wildtype"))

lapply(averageRowIDs, )


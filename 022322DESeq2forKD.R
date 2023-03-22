ReportersSubJGIAllcountsKOoI2 <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReportersSubJGIAllcountsKOoI2.csv", row.names = 1, check.names = FALSE, header=TRUE, sep=",")
ReportersSubJGIAllcountsKOoImatrix2 <- data.matrix(ReportersSubJGIAllcountsKOoI2)
ReportersSubJGIAllcountsKOoIcountsOnly2plus <- data.matrix(ReportersSubJGIAllcountsKOoImatrix2[ ,6:ncol(ReportersSubJGIAllcountsKOoI2)])+0.00001

ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly2plus, lengths = ReportersSubJGIAllcountsKOoI2[,5])

metadata <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/metadata_ReportersSubJGIAllCountsKOoI.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
rownames(metadata) <- make.names(metadata[,1], unique = TRUE)
colnames(ReportersSubJGIAllcountsKOoIcountsOnly2plus) <- rownames(metadata)
ReportersSubJGIAllcountsKOoIcountsOnly2plus  <- data.matrix(ReportersSubJGIAllcountsKOoIcountsOnly2plus)

ReportersnonClusteredKOoI<- pheatmap(Averaged_Orderd_KO_data[,-29], color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(200), cellwidth = 20, cellheight = 30, legend=T, show_rownames=T, show_colnames=T, fontsize_col=15, fontsize_row = 10, treeheight_row=0, treeheight_col=20, height = 1.5, width = 1.5)


ddsAVG <- DESeqDataSetFromMatrix(countData = round(ReportersSubJGIAllcountsKOoIcountsOnly2plus), colData = metadata, design = ~condition)
ddsAVG$condition <- relevel(ddsAVG$condition, ref = "wildtype")
ddsAVG$condition <- factor(ddsAVG$condition)
ddsAVG$type <- factor(ddsAVG$type)
ddsAVG1 <- DESeq(ddsAVG)
summary(ddsAVG1)

averageRowIDs =c("set7", "isw", "suz12", "eed", "crf4.3", "ncu00548", "ncu06787", "ncu06788", "H3K4R", "ncu06788", "ncu04017", "ncu00548", "cac3", "cac2", "ncu00423", "ncu00423", "ncu00548", "set2", "ncu09120", "ncu03481", "ncu03481", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017", "rtt109", "suz12", "wildtype")

averageRowIDs_no_wt =c("set7", "isw", "suz12", "eed", "crf4.3", "ncu00548", "ncu06787", "ncu06788", "H3K4R", "ncu06788", "ncu04017", "ncu00548", "cac3", "cac2", "ncu00423", "ncu00423", "ncu00548", "set2", "ncu09120", "ncu03481", "ncu03481", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017", "rtt109", "suz12")

results(ddsAVG1, contrast=c("condition", "set7", "wildtype"))

deseq2_results_writer <- function(strain_name, deseq2_obj, output_dir){
  
  mutant_result_1 <- results(deseq2_obj, alpha=.01, contrast=c("condition", strain_name, "wildtype"))
  
  return(mutant_result_1)
}


x <- lapply(averageRowIDs_no_wt, deseq2_results_writer, ddsAVG1, "test")
str(x)

csv <- function(log2foldchange, x, output_column())
  write.csv(x, file = "testing.csv")

l2fc <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/022123REPORTERSsubJGIallcountsKOoI/log2fc.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE, sep=",")
l2fc <- as.matrix(l2fc)
barp <- barplot(l2fc[6,], main="NCU09953", las = 3,
                ylab="Log 2 Fold Change", col = brewer.pal(12, name="Set3"))
barp <- barplot(l2fc[4,], main="NCU07149", las = 3,
                ylab="Log 2 Fold Change", col = brewer.pal(12, name="Set3"))
barp <- barplot(abs(l2fc[3,]), main="NCU06889", las = 3,
                ylab="Log 2 Fold Change", col = brewer.pal(12, name="Set3"))
barp <- barplot(l2fc[5,], main="NCU08085", las = 3,
                ylab="Log 2 Fold Change", col = brewer.pal(12, name="Set3"))


ddsAVG1name <- write.table(resultsNames(ddsAVG1))

#calculateTPM
ReportersSubJGIAllcountsKOoI2TPM <- calculateTPM(ReportersSubJGIAllcountsKOoIcountsOnly2plus, lengths = ReportersSubJGIAllcountsKOoIcountsOnly2plus[,5])
##convert to matrix again
ReportersSubJGIAllcountsKOoI2TPM  <- data.matrix(ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM)
Ordered_KO_data <- cbind(ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,1:3],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,4:5],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,6:8],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,9:11],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,12:14],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,15:17],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,18:20],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,21:23],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,24:26],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,27:29],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,30:32],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,33:35],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,36:38],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,39:41],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,42:44],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,45:47],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,48:50],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,51:53],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,54:56],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,57:58],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,59:61],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,62:64],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,65:67],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,68:70],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,71:73],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,74:76],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,77:79],
                         ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,80:82],ReportersSubJGIAllcountsKOoIcountsOnly2plusTPM[,83:85])
Averaged_Orderd_KO_data <- cbind(rowMeans(Ordered_KO_data[,1:3], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,4:5], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,6:8], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,9:11], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,12:14], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,15:17], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,18:20], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,21:23], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,24:26], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,27:29], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,30:32], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,33:35], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,36:38], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,39:41], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,42:44], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,45:47], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,48:50], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,51:53], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,54:56], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,57:58], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,59:61], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,62:64], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,65:67], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,68:70], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,71:72], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,74:76], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,77:79], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,80:82], na.rm = TRUE),
                                 rowMeans(Ordered_KO_data[,83:85], na.rm = TRUE))
averageRowIDs =c("set7", "isw", "suz12a", "eed", "crf4.3", "ncu00548a", "ncu06787", "ncu06788a", "H3K4R", "ncu06788A", "ncu04017a", "ncu00548A", "cac3", "cac2", "ncu00423a", "ncu00423A", "ncu00548a", "set2", "ncu09120", "ncu03481a", "ncu03481A", "ncu02695", "ncu03461", "ncu03461", "nst1", "ncu04017A", "rtt109", "suz12A", "wildtype")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

barp <- barplot(Averaged_Orderd_KO_data[3,-29], main="NCU06889", las = 3,
                ylab="Average TPM", col = brewer.pal(12, name="Set3"), ylim=c(0, 2e+05))
write.csv(Averaged_Orderd_KO_data, file="average_TPMperKO.csv")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
library("scater")
library("dplyr")
library("tidyverse")
library("DESeq2")
library("ggplot2")
library("scales")
library(readr)
library("pheatmap")
library("RColorBrewer")

Allcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JgiAllSampleCounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
K27SilentAcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JGI_ListofK27SilentGenesNCUs.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
SRRtoGENE <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/RNASeqAllJGI_KOonly.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
#selected <- read.table("/Users/ry00555/Desktop/SelectedKO.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
#selected2 <- read.table("/Users/ry00555/Desktop/SelectedKO2.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
#selectedKOmetadata2 <- read.table("/Users/ry00555/Desktop/SelectedKO2metadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
#selectedKOmetadata <- read.table("/Users/ry00555/Desktop/SelectedKOmetadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
Allcountsmatrix <- data.matrix(Allcounts)
Allcounts_countsOnly <- data.matrix(Allcountsmatrix[ ,6:ncol(Allcountsmatrix)])
colnames(Allcounts_countsOnly) <- SRRtoGENE[,1]
Allcounts_TPM <- calculateTPM(Allcounts_countsOnly, lengths = Allcountsmatrix[,5])
allDataTPM  <- data.matrix(Allcounts_TPM)
PRC2TargetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%K27SilentAcounts[,1])


Topselected2 <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Top20KOwPRC2GeneChanges.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
selected_cols2 <- allDataTPM[,Topselected2$V1]
GenesWithChanges_SelectedKO2 <- subset(selected_cols2, (rowSums(selected_cols2) > 0))


selected2_counts <- Allcounts_countsOnly[,Topselected2$V1]
cts_select2 <- as.matrix(selected2_counts)
cts_select2new_numeric <- matrix(as.numeric(cts_select2), ncol = 23)
dimnames(cts_select2new_numeric) <- list(rownames(cts_select2), colnames(cts_select2))
Top20_shortened <- as.data.frame(as.matrix(cbind(cts_select2new_numeric[,1:3], cts_select2new_numeric[,4:5], cts_select2new_numeric[,6:8], cts_select2new_numeric[,9], cts_select2new_numeric[,10:12], cts_select2new_numeric[,13:15], cts_select2new_numeric[,16:18], cts_select2new_numeric[,19:21], cts_select2new_numeric[,22], cts_select2new_numeric[,23])))

TopselectedKOmetadata2 <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Top20KOwPRC2GeneChangesMetadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
rownames(TopselectedKOmetadata2) <- make.names(TopselectedKOmetadata2[,1], unique = TRUE)
cts_select2 <- Top20_shortened
colnames(cts_select2) <- rownames(TopselectedKOmetadata2)

TopselectedKOmetadata2$condition <- factor(TopselectedKOmetadata2$condition)
TopselectedKOmetadata2$type <- factor(TopselectedKOmetadata2$type)
colnames(cts_select2) <- rownames(TopselectedKOmetadata2)
all(rownames(TopselectedKOmetadata2) %in% colnames(cts_select2))
all(rownames(TopselectedKOmetadata2) == colnames(cts_select2))

dds3<- DESeqDataSetFromMatrix(countData = cts_select2,
                              colData = TopselectedKOmetadata2,
                              design = ~ condition)
dds3

keep <- rowSums(counts(dds3))>= 10
dds3<- dds3[keep,]
dds3$condition <- relevel(dds3$condition, ref = "wildtype")

dds3 <- DESeq(dds3)
plotDispEsts(dds3)
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/VolcanoPlotTop20/Sig")
ggsave(filename = "DeSeq_TopelectKO2_DispersionEstimate.pdf", plot = dds3, dpi=4200)

alpha = 0.05
padj.cutoff <- 0.05
res_isw <- results(dds3, alpha=alpha, contrast=c("condition", "isw", "wildtype"))
res_ncu06787 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu06787", "wildtype"))
res_rtt109 <- results(dds3, alpha=alpha, contrast=c("condition", "rtt109", "wildtype"))
res_cac3 <- results(dds3, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_H3K4R <- results(dds3, alpha=alpha, contrast=c("condition", "H3K4R", "wildtype"))
res_ncu06788 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu06788", "wildtype"))
res_ncu04117 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu04017", "wildtype"))
res_ncu00548 <- results(dds3, alpha=alpha, contrast=c("condition", "ncu00548", "wildtype"))

res_sig_isw <- subset(res_isw, padj < 0.05)
res_sig_ncu06787 <- subset(res_ncu06787, padj < 0.05)
res_sig_rtt109 <- subset(res_rtt109, padj < 0.05)
res_sig_cac3 <- subset(res_cac3, padj < 0.05)
res_sig_H3K4R <- subset(res_H3K4R, padj < 0.05)
res_sig_ncu06788 <- subset(res_ncu06788, padj < 0.05)
res_sig_ncu041017 <- subset(res_ncu04117, padj < 0.05)
res_sig_ncu00548 <- subset(res_ncu00548, padj < 0.05)

write.csv(as.data.frame(res_sig_isw), file = "res_sig_isw.csv")
write.csv(as.data.frame(res_sig_ncu06787), file = "res_sig_ncu06787.csv")
write.csv(as.data.frame(res_sig_rtt109), file = "res_sig_rtt109.csv")
write.csv(as.data.frame(res_sig_cac3), file = "res_sig_cac3.csv")
write.csv(as.data.frame(res_sig_H3K4R), file = "res_sig_H3K4R.csv")
write.csv(as.data.frame(res_sig_ncu06788), file = "res_sig_ncu06788.csv")
write.csv(as.data.frame(res_sig_ncu041017), file = "res_sig_ncu04017.csv")
write.csv(as.data.frame(res_sig_ncu00548), file = "res_sig_ncu00548.csv")


list_of_sigfiles <- list.files(path = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/VolcanoPlotTop20/Sig ",
                               recursive = TRUE,
                               pattern = ".csv$")

induced_genes_sig2 <- readr::read_csv(list_of_sigfiles, id = "file_name")
#l2fc_selectedKO2_sigdata <- data.frame(pivot_wider(data = induced_genes_sig2, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))

l2fc_Top20 <- data.frame(pivot_wider(data = induced_genes_sig2, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))

rownames(l2fc_Top20) <- l2fc_Top20[,1]
l2fc_Top20 <- l2fc_Top20[,-1]
l2fc_Top20[is.na(l2fc_Top20)]<-0
GenesWithChanges2 <- subset(l2fc_Top20, (rowSums(l2fc_Top20) > 0))
l2fc_Top20_PRC2Targetonly <- subset(l2fc_Top20, rownames(l2fc_Top20)%in%K27SilentAcounts[,1])
l2fc_Top20_PRC2Targetonly <- subset(l2fc_Top20_PRC2Targetonly, (rowSums(l2fc_Top20_PRC2Targetonly) > 0))


breaks <- c(seq(0, 15, by=0.1))
V2Top20_PRC2Target_Heatmap<- pheatmap(l2fc_Top20_PRC2Targetonly, breaks=breaks, color = colorRampPalette((brewer.pal(n = 9, name="YlOrRd")))(200), cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)
V2Top20_PRC2Target_Heatmaps_plot<-V2Top20_PRC2Target_Heatmap[[4]]
ggsave(filename = "Top20_PRC2Target_HeatmapNewColorSameGenes.pdf", plot = V2Top20_PRC2Target_Heatmaps_plot, dpi=4200)


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
induced_genes_isw <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/VolcanoPlotTop20/Sig /Iswi/res_sig_isw.csv")
rownames(induced_genes_isw) <- induced_genes_isw[,1]
induced_genes_isw <- induced_genes_isw[,-1]
IswGenesWithChanges2 <- subset(induced_genes_isw, (rowSums(induced_genes_isw) > 0))
l2fc_isw_PRC2Targetonly <- subset(IswGenesWithChanges2, rownames(IswGenesWithChanges2)%in%K27SilentAcounts[,1])


EnhancedVolcano(l2fc_isw_PRC2Targetonly,
                lab = rownames(l2fc_isw_PRC2Targetonly),
                title = 'isw versus wildtype',
                pCutoff = 10e-32,
                x = 'log2FoldChange',
                y = 'pvalue')


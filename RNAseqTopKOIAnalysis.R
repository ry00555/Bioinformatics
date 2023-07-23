setwd('/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts')
library(rtables)
library("scater")
library("dplyr")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("ggplot2")
library("scales")
library(readr)


AllcountsKOname <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JgiAllSampleCountswKOname.csv", header=TRUE, stringsAsFactors=TRUE, row.names=1, check.names=FALSE, sep=",")
AllcountsKOnameMATRIX <- data.matrix(AllcountsKOname)
AllcountsKOnameCountsOnly <- (AllcountsKOnameMATRIX[ ,6:ncol(AllcountsKOnameMATRIX)])
AllcountsKOnameCountsOnly <- data.matrix(AllcountsKOnameMATRIX[ ,6:ncol(AllcountsKOnameMATRIX)])

K27SilentAcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JGI_ListofK27SilentGenesNCUs.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
ReporterGenes <- read.table('/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/ReporterGeneNames.txt', header=FALSE, check.names=FALSE, sep="\t")

calculateTPM
Allcounts_TPM <- calculateTPM(AllcountsKOnameCountsOnly, lengths = AllcountsKOnameMATRIX[,5])
Allcounts_TPM  <- data.matrix(Allcounts_TPM)

PRC2TargetTPM <- subset(Allcounts_TPM, rownames(Allcounts_TPM)%in%K27SilentAcounts[,1])
ReportersTPM <- subset(Allcounts_TPM, rownames(Allcounts_TPM)%in%ReporterGenes[,1])
write.csv(PRC2TargetTPM, file="PRC2TargetTPM.csv")

PRC2TargetTPMheatmap<- pheatmap(PRC2TargetTPM, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 50, width = 50)
ReportersTPMheatmap<- pheatmap(ReportersTPM, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = 20,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=10, height = 80, width = 50)

Top50 <-PRC2TargetTPMheatmap$tree_col$order[1:50]
Top50names <- PRC2TargetTPMheatmap$tree_col$labels[Top50]
PRc2TargetTop50KO <- PRC2TargetTPM[,Top50names]
ImportPRc2TargetTop50KO  <- data.matrix(PRc2TargetTop50KO)
ImportPRc2TargetTop50KOmorethan0 <- subset(ImportPRc2TargetTop50KO, (rowSums(ImportPRc2TargetTop50KO) > 0))
write.csv(ImportPRc2TargetTop50KOmorethan0, file="PRC2TargetTop50TPM.csv")
heatmap<- pheatmap(ImportPRc2TargetTop50KOmorethan0, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100), breaks=breaks, cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 5, width = 2.5)

ReportersTop50 <-ReportersTPMheatmap$tree_col$order[1:50]
ReportersTop50names <- ReportersTPMheatmap$tree_col$labels[Top50]
ReportersTargetTop50KO <- ReportersTPM[,ReportersTop50names]
ReportersTargetTop50KO  <- data.matrix(ReportersTargetTop50KO)
write.csv(ReportersTargetTop50KO, file="ReportersTargetTop50KO.csv")
breaks <- c(seq(0, 4, by=0.05))
heatmap<- pheatmap(ReportersTargetTop50KO, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100), breaks=breaks, cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, lablegend='Average TPM', show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=5, treeheight_col=10, height = 5, width = 2.5)

KONames <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/KONames.txt", header=TRUE, stringsAsFactors=TRUE, check.names=FALSE, sep="\t")          
KOoICounts <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/KOoIAllCounts.csv", header=TRUE,  stringsAsFactors=FALSE, check.names=FALSE, sep=",", row.names=1)
KOoICountsMatrix <- data.matrix(KOoICounts)
KOoICountsOnlyMatrix <- data.matrix(KOoICountsMatrix[ ,6:ncol(KOoICountsMatrix)])

KOoITPM <- data.matrix(calculateTPM(KOoICountsOnlyMatrix, lengths = KOoICountsMatrix[,5]))
KOoITPMFilterNegative <- subset(KOoITPM, (rowSums(KOoITPM) > 0))
KOoITPMheatmap<- pheatmap(KOoITPMFilterNegative, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 50, width = 50)

KOoIPRC2TargetTPM <- subset(KOoITPMFilterNegative, rownames(KOoITPMFilterNegative)%in%K27SilentAcounts[,1])
KOoIPRC2TargetTPMheatmap<- pheatmap(KOoIPRC2TargetTPM, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 50, width = 50)

KOoIReportersTPM <- subset(KOoITPMFilterNegative, rownames(KOoITPMFilterNegative)%in%ReporterGenes[,1])
KOoIReportersTPMheatmap<- pheatmap(KOoIReportersTPM, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 50, width = 50)

KOoIMeta <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/PRC2TargetKOoI.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
rownames(KOoIMeta) <- make.names(KOoIMeta[,1], unique = TRUE)
colnames(KOoICountsOnlyMatrix) <- rownames(KOoIMeta)
colnames(KOoICountsOnlyMatrix) == rownames(KOoIMeta)

KOoICountsOnlyMatrix_numeric <- matrix(as.numeric(KOoICountsOnlyMatrix), ncol = 103)
colnames(KOoICountsOnlyMatrix_numeric) <- rownames(KOoIMeta)
rownames(KOoICountsOnlyMatrix_numeric) <- rownames(KOoICountsOnlyMatrix)
colnames(KOoICountsOnlyMatrix_numeric) == colnames(KOoICountsOnlyMatrix)

dimnames(KOoICountsOnlyMatrix_numeric) <- list(rownames(KOoICountsOnlyMatrix), colnames(KOoICountsOnlyMatrix))
KOoIShortened <- as.data.frame(as.matrix(cbind(KOoICountsOnlyMatrix_numeric[,1:3], KOoICountsOnlyMatrix_numeric[,4:6], KOoICountsOnlyMatrix_numeric[,7:9], KOoICountsOnlyMatrix_numeric[,10:12], KOoICountsOnlyMatrix_numeric[,13:15], KOoICountsOnlyMatrix_numeric[,16:18], KOoICountsOnlyMatrix_numeric[,19:21], KOoICountsOnlyMatrix_numeric[,22:24], KOoICountsOnlyMatrix_numeric[,25:27], KOoICountsOnlyMatrix_numeric[,28:30], 
                                                 KOoICountsOnlyMatrix_numeric[,31:33], KOoICountsOnlyMatrix_numeric[,34:36], KOoICountsOnlyMatrix_numeric[,37:39], KOoICountsOnlyMatrix_numeric[,40:42], KOoICountsOnlyMatrix_numeric[,43:45], KOoICountsOnlyMatrix_numeric[,46:47], KOoICountsOnlyMatrix_numeric[,48:49], KOoICountsOnlyMatrix_numeric[,50:52], KOoICountsOnlyMatrix_numeric[,53:54], KOoICountsOnlyMatrix_numeric[,55:59], 
                                                 KOoICountsOnlyMatrix_numeric[,60:62], KOoICountsOnlyMatrix_numeric[,63:65], KOoICountsOnlyMatrix_numeric[,66:68], KOoICountsOnlyMatrix_numeric[,69:71], KOoICountsOnlyMatrix_numeric[,72:74], KOoICountsOnlyMatrix_numeric[,75:77], KOoICountsOnlyMatrix_numeric[,78:80], KOoICountsOnlyMatrix_numeric[,81:83], KOoICountsOnlyMatrix_numeric[,84:86], KOoICountsOnlyMatrix_numeric[,87:89], 
                                                 KOoICountsOnlyMatrix_numeric[,90:92], KOoICountsOnlyMatrix_numeric[,93:95], KOoICountsOnlyMatrix_numeric[,96:98], KOoICountsOnlyMatrix_numeric[,99:103] )))
KOoIMetaRowIDs =c("0.17%arginine", "NCU00202A", "NCU00269 NCU07496 NCO00119", "NCU00423A", "NCU00548a", "NCU01238a", "NCU01238A",
              "NCU01635K36R", "NCU01635K4R","NCU01888", "NCU02017a", "NCU02695", "NCU03060", "NCU03073", "NCU03481a",
              "NCU03481A", "NCU03875", "NCU04017", "NCU04198A", "NCU04737","NCU05300A",  "NCU05460a", "NCU05460A",
              "NCU06679", "NCU06787A", "NCU06788a", "NCU06788A", "NCU07496", "NCU08357", "NCU09120",
              'NCU09825A',"Wildtype")

KOoIMetaRowIDsnoWT =c("0.17%arginine", "NCU00202A", "NCU00269 NCU07496 NCO00119", "NCU00423A", "NCU00548a", "NCU01238a", "NCU01238A",
                      "NCU01635K36R", "NCU01635K4R","NCU01888", "NCU02017a", "NCU02695", "NCU03060", "NCU03073", "NCU03481a",
                      "NCU03481A", "NCU03875", "NCU04017", "NCU04198A", "NCU04737","NCU05300A",  "NCU05460a", "NCU05460A",
                      "NCU06679", "NCU06787A", "NCU06788a", "NCU06788A", "NCU07496", "NCU08357", "NCU09120",
                      'NCU09825A')



KOoIdds2 <- DESeqDataSetFromMatrix(countData = KOoIShortened,
                                   colData = KOoIMeta,
                                   design = ~ Type)

KOoIMeta$Condition <- factor(KOoIMeta$Condition)
KOoIMeta$Type <- factor(KOoIMeta$Type)

KOoIdds2$Type <- relevel(KOoIdds2$Type, ref = "Wildtype")
KOoIdds2 <- DESeq(KOoIdds2)
plotDispEsts(KOoIdds2)
resultsNames(KOoIdds2)
set7res <- write.table(results(KOoIdds2, contrast=c("Type", "NCU07496", "Wildtype")))

deseq2_results_writer <- function(strain_name, KOoIdds2, output_dir){
  
  mutant_result_1 <- results(KOoIdds2, alpha=.005, contrast=c("Type", strain_name, "Wildtype"))
  
  return(mutant_result_1)
}


x <- lapply(KOoIMetaRowIDsnoWT, deseq2_results_writer, KOoIdds2, "test")
str(x)

KOoIDeSeq2L2FC <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/KOoIL2FC_0306023.csv", header = TRUE, stringsAsFactors = FALSE, sep=",", row.names=1)
KOoIDeSeq2pValue <- read.csv("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/KOoIpValue_030623.csv", header = TRUE, stringsAsFactors = FALSE, sep=",", row.names=1)
rownames(KOoIDeSeq2L2FC) == rownames(KOoIDeSeq2pValue)
KOoIDeSeq2L2FCmatrix <- data.matrix(KOoIDeSeq2L2FC)
KOoIDeSeq2pValuematrix <- data.matrix(KOoIDeSeq2pValue)


#Attempting to plot L2FC against p-value for all strains with each dot representing a gene. 
df <-data.frame(x=KOoIDeSeq2pValuematrix, y=KOoIDeSeq2L2FCmatrix)
dat <- data.frame(a=as.vector(KOoIDeSeq2pValuematrix[upper.tri(KOoIDeSeq2pValuematrix)]),
                  b=as.vector(KOoIDeSeq2L2FCmatrix[upper.tri(KOoIDeSeq2L2FCmatrix)]))
plot(dat$a,dat$b)

GenesfromDeSeq <- rownames(KOoIDeSeq2L2FC)
print(GenesfromDeSeq)

KOoIDeSeq2L2FCPRC2Target <- data.matrix(subset(KOoIDeSeq2L2FC, rownames(KOoIDeSeq2L2FC)%in%K27SilentAcounts[,1]))
#Check if there are any NA values that need to be changed to 0 before plotting 
any(is.na(KOoIDeSeq2L2FCPRC2Target))
which(is.na(KOoIDeSeq2L2FCPRC2Target))
KOoIDeSeq2L2FCPRC2Target[is.na(KOoIDeSeq2L2FCPRC2Target)] <- 0.1
#No rows with all NAs or zero variance like you said, but if you do the dist calculation, there are NAs in some of the entries, indicating between some rows, it's not possible to calculate euclidean distances. You need to the euclidean distance matrix to have no NAs to do clustering:

sum(is.na(as.matrix(dist(KOoIDeSeq2L2FCPRC2Target))))
pheatmap(KOoIDeSeq2L2FCPRC2Target,  width = 50, height = 50, show_rownames=T, show_colnames=T, clustering_method="centroid", cellheight = 6, fontsize_row = 5)
KOoIDeSeq2L2FCPRC2TargetHeatmap<- pheatmap(KOoIDeSeq2L2FCPRC2Target2, color = colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=T, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 5, width = 2.5)
plot(KOoIDeSeq2L2FCPRC2Target)

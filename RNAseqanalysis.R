for i in SRR8269612

do 
wget -O $i.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$i&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,sample_alias,sample_title&format=tsv&download=true&limit=0"
done

setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq")

#Calculate TPM
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
```{r calculate TPM}
library("scater")
library("dplyr")
library("tidyverse")

Allcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JgiAllSampleCounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
Allcountsmatrix <- data.matrix(Allcounts)
Allcounts_countsOnly <- data.matrix(Allcountsmatrix[ ,6:ncol(Allcountsmatrix)])
K27SilentAcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JGI_ListofK27SilentGenesNCUs.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
SRRtoGENE <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/RNASeqAllJGI_KOonly.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
colnames(Allcounts_countsOnly) <- SRRtoGENE[,1]
iswset7 <- read.table("/Users/ry00555/Desktop/iswset7.txt", header=FALSE)


calculateTPM
Allcounts_TPM <- calculateTPM(Allcounts_countsOnly, lengths = Allcountsmatrix[,5])
##convert to matrix again
allDataTPM  <- data.matrix(Allcounts_TPM)

#read in PRC2 target genes table
K27SilentAcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI_ListofK27SilentGenesNCUs.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
PRC2TargetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%K27SilentAcounts[,1])
iswset7only <- select(c(PRC2TargetTPM, 'crf4-1', 'crf4-1-3'))

iswset7only <- subset(PRC2TargetTPM, rownames(PRC2TargetTPM)%in%iswset7[,1])

#Zack's commands 
countsonly<- as.matrix(raw_counts_df[,6:ncol(raw_counts_df)])
as.numeric(raw_counts_df$Length)
all_tpm <- calculateTPM(countsonly, effective_length = raw_counts_df$Length)

Prc2targetTPM <- subset(all_tpm, rownames(all_tpm)%in%Prc2targets[,1])



#rename columns with human readable sample names

###THIS NEXT SECTION REQUIRED MANUAL INPUT INTO A SPREADSHEET OR TEXT EDITOR
#make a list of just the gene names 

SRRtoGENE <- read.table("/Users/ry00555/Desktop/RNASeqAllJGI_KOonly.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
colnames(allDataTPM) <- SRRtoGENE[,1]
colnames(PRC2TargetTPM) <- SRRtoGENE[,1]
colnames(Allcounts_countsOnly) <- SRRtoGENE[,1]
genelist <- rownames(allDataTPM)
#make a list of just the SRRIDs
SRRnames <-colnames(allDataTPM)

#make a table of just SRR names
write.table(SRRnames, file="SRRnames.txt")

##Keep Only Genes that are expressed in at least one sample
#GenesWithChanges <- subset(Prc2targetTPM, (rowSums(Prc2targetTPM) > 0))
GenesWithChanges <- subset(PRC2TargetTPM, (rowSums(PRC2TargetTPM) > 0))
PRC2TargetTPM <- subset(GenesWithChanges, colnames(GenesWithChanges)%in%K27SilentAcounts[,1])

```


#Make a heatmap

```{r make a heatmap}

#make a heatmap parameters
library("pheatmap")
library("RColorBrewer")

heatmap<- pheatmap(GenesWithChanges, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot<-heatmap[[4]]

#save a heatmap file
ggsave(filename = "./allStrains_GenesWChanges.pdf", plot = heatmap_plot, dpi=600)

dev.off()




```

#Make a smaller heatmap with the 100 samples with the highest expression levels of PRC2 induced genes

```{r}
##gen selected names from heatmap
samplenumbers <-heatmap$tree_col$order[1:20]
samplenames <- heatmap$tree_col$labels[samplenumbers]

##Keep Only 100 Samples with highest expression levels
#Prc2targetTPM <- subset(all_tpm, rownames(all_tpm)%in%Prc2targets[,1])
SelectedSamplesTop20 <- PRC2TargetTPM[,samplenames]
GenesWithChanges_selectedSamples <- subset(SelectedSamplesTop20, (rowSums(SelectedSamplesTop20) > 0))

write.table(samplenames, file="./SampleName1to50.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

heatmap2Top20genename<- pheatmap(GenesWithChanges_selectedSamples, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 10)
heatmap2Top20genename<-heatmap2Top20genename[[4]]
ggsave(filename = "heatmapTop20_GenesWChanges_PRC2Target.pdf", plot = heatmap2Top20genename, dpi=4200)


GenesWithChanges_selectedSamples <- subset(allDataTPM, colnames(allDataTPM)%in%selected$V1)
selected <- subset(PRC2TargetTPM, colnames(PRC2TargetTPM) %in% selected2)

heatmapTop20_GenesWChanges_PRC2Target <- pheatmap(GenesWithChanges_selectedSamples, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = T, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 10)
heatmapTop20_GenesWChanges_PRC2Target<-heatmapTop20_GenesWChanges_PRC2Target[[4]]
ggsave(filename = "heatmapTop20_GenesWChanges_PRC2Target.pdf", plot = heatmapTop20_GenesWChanges_PRC2Target, dpi=4200)
#select KO strains of interest 
selected <- read.table("/Users/ry00555/Desktop/SelectedKO.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
selected_cols2 <- allDataTPM[,selected2$V1]
GenesWithChanges <- subset(selected_cols2, (rowSums(selected_cols2) > 0))
heatmap3<- pheatmap(GenesWithChanges, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, scale="row", cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)
heatmap_plot3<-heatmap3[[4]]
ggsave(filename = "Heatmap_SelectedKO2.pdf", plot = heatmap_plot3, dpi=4200)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot2<-heatmap2[[4]]

#save a heatmap file
ggsave(filename = "./Heatmap_sampels1to20.pdf", plot = heatmap_plot2, dpi=4200)


dev.off()


```


#Make a boxplot

```{r make a boxplot}

##make melted dataframe for boxplots
library(reshape2)
meltedPRC2targetData <- melt(Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))
```````````
#DEseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("scales")
BiocManager::install("IHW")
library("DESeq2")
library("ggplot2")
library("scales")
library(readr)

## Importing featurecounts matrix
Read in the featurecounts output, convert to a matrix, and remove unneccesary columns to include ONLY the count numbers, then convert columns to numeric.
```{r cts}
selected_counts <- Allcounts_countsOnly[,selected$V1]
cts <- as.matrix(selected_counts)
cts_numeric <- matrix(as.numeric(cts), ncol = 15)
dimnames(cts_numeric) <- list(rownames(cts), colnames(cts))

cts_shortened <- as.data.frame(as.matrix(cbind(cts_numeric[,1:3], cts_numeric[,4:6], cts_numeric[,7:9], cts_numeric[,10:12], cts_numeric[,13:15])))

#import metadata for select KO strains 
selectedKOmetadata <- read.table("/Users/ry00555/Desktop/SelectedKOmetadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
rownames(selectedKOmetadata) <- make.names(selectedKOmetadata[,1], unique = TRUE)
cts2 <- cts_shortened
colnames(cts2) <- rownames(selectedKOmetadata)


selectedKOmetadata$condition <- factor(selectedKOmetadata$condition)
selectedKOmetadata$type <- factor(selectedKOmetadata$type)

all(rownames(selectedKOmetadata) %in% colnames(cts2))
all(rownames(selectedKOmetadata) == colnames(cts2))


##create DEseq dataset
dds2<- DESeqDataSetFromMatrix(countData = cts2,
                              colData = selectedKOmetadata,
                              design = ~ condition)
dds2

##pre-filteringreads <10
keep <- rowSums(counts(dds2))>= 10
dds2<- dds2[keep,]

#specify level of comparison
dds2$condition <- relevel(dds2$condition, ref = "FGSC4200")


#Now, run DEseq and plot dispesion estimates to see the quality of the fit. 

```{r coldata}
dds2 <- DESeq(dds2)
plotDispEsts(dds2)
```
Export genes with statistically significant differential expression when compared against the wild type.

```{r coldata}
alpha = 0.05

res_S540 <- results(dds2, alpha=alpha, contrast=c("condition", "S540" , "FGSC4200"))
res_S238 <- results(dds2, alpha=alpha, contrast=c("condition", "S238" , "FGSC4200"))
res_FGSC12770 <- results(dds2, alpha=alpha, contrast=c("condition", "FGSC12770" , "FGSC4200"))
res_FGSC14852 <- results(dds2, alpha=alpha, contrast=c("condition", "FGSC14852" , "FGSC4200"))


write.csv(as.data.frame(res_S540), file = "res_S540.csv")
write.csv(as.data.frame(res_S238), file = "res_S238.csv")
write.csv(as.data.frame(res_FGSC12770), file = "res_FGSC12770.csv")
write.csv(as.data.frame(res_FGSC14852), file = "res_FGSC14852.csv")


Subset to only export significant differentially expressed genes

```{r}
res_sig_S540 <- subset(res_S540, padj < 0.05)
res_sig_S238 <- subset(res_S238, padj < 0.05)
res_sig_FGSC12770 <- subset(res_FGSC12770, padj < 0.05)
res_sig_FGSC14852 <- subset(res_FGSC14852, padj < 0.05)


write.csv(as.data.frame(res_sig_S540), file = "res_sig_S540.csv")
write.csv(as.data.frame(res_sig_S238), file = "res_sig_S238.csv")
write.csv(as.data.frame(res_sig_FGSC12770), file = "res_sig_FGSC12770.csv")
write.csv(as.data.frame(res_sig_FGSC14852), file = "res_sig_FGSC14852.csv")

setwd("/Users/ry00555/Desktop/SelectedKO11122")
list_of_files <- list.files(path = "/Users/ry00555/Desktop/SelectedKO11122",
                            recursive = TRUE,
                            pattern = ".csv$")

induced_genes <- readr::read_csv(list_of_files, id = "file_name")


##convert to wide format
l2fc_data <- data.frame(pivot_wider(data = induced_genes, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_data) <- l2fc_data[,1]
l2fc_data <- l2fc_data[,-1]
PRC2Targetgenesl2fc <-subset(l2fc_data, rownames(l2fc_data)%in%K27SilentAcounts[,1])
PRC2Targetgenesl2fc <- PRC2Targetgenesl2fc[,-(5:8)]

#heat map with log 2 fold change 
selected <- read.table("/Users/ry00555/Desktop/SelectedKO.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
selected_cols <- allDataTPM[,selected$V1]
GenesWithChanges <- subset(selected_cols, (rowSums(selected_cols) > 0))

heatmap4<- pheatmap(l2fc_data, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)
heatmap_plot4<-heatmap4[[4]]

heatmap5<- pheatmap(PRC2Targetgenesl2fc, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)

ggsave(filename = "Heatmap_Selected.pdf", plot = heatmap_plot3, dpi=4200)

``````````````````

```````````
ON 11/1/2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("scales")
BiocManager::install("IHW")
library("DESeq2")
library("ggplot2")
library("scales")
library(readr)

#for selected KO mutants 2 
selected2 <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/SignificantSelectedKO11122/SelectedKO2.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
selected_cols2 <- allDataTPM[,selected2$V1]
GenesWithChanges_SelectedKO2 <- subset(selected_cols2, (rowSums(selected_cols2) > 0))


## Importing featurecounts matrix
Read in the featurecounts output, convert to a matrix, and remove unneccesary columns to include ONLY the count numbers, then convert columns to numeric.
```{r cts}
selected2_counts <- Allcounts_countsOnly[,selected2$V1]
cts_select2 <- as.matrix(selected2_counts)
cts_select2new_numeric <- matrix(as.numeric(cts_select2), ncol = 10)
dimnames(cts_select2new_numeric) <- list(rownames(cts_select2), colnames(cts_select2))
cts_select2_newshortened <- as.data.frame(as.matrix(cbind(cts_select2new_numeric[,1:2], cts_select2new_numeric[,3:4], cts_select2new_numeric[,5:6], cts_select2new_numeric[,7:8], cts_select2new_numeric[,9:10])))

#import metadata for select KO strains 
selectedKOmetadata2 <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/SignificantSelectedKO11122/SelectedKO2metadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
rownames(selectedKOmetadata2) <- make.names(selectedKOmetadata2[,1], unique = TRUE)
cts_select2 <- cts_select2_newshortened
colnames(cts_select2) <- rownames(selectedKOmetadata2)

````````````
FOR REFERENCE ONLY
selected_counts <- Allcounts_countsOnly[,selected$V1]
cts <- as.matrix(selected_counts)
cts_numeric <- matrix(as.numeric(cts), ncol = 15)
dimnames(cts_numeric) <- list(rownames(cts), colnames(cts))

cts_shortened <- as.data.frame(as.matrix(cbind(cts_numeric[,1:3], cts_numeric[,4:6], cts_numeric[,7:9], cts_numeric[,10:12], cts_numeric[,13:15])))

#import metadata for select KO strains 
selectedKOmetadata <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/SignificantSelectedKO11122/SelectedKO2metadata.txt", header=TRUE, row.names =1, stringsAsFactors=FALSE, sep="\t")
rownames(selectedKOmetadata) <- make.names(selectedKOmetadata[,1], unique = TRUE)
cts2 <- cts_shortened
colnames(cts2) <- rownames(selectedKOmetadata)
````````````````````````````


selectedKOmetadata2$condition <- factor(selectedKOmetadata2$condition)
selectedKOmetadata2$type <- factor(selectedKOmetadata2$type)
colnames(cts_select2) <- rownames(selectedKOmetadata2)
all(rownames(selectedKOmetadata2) %in% colnames(cts_select2))
all(rownames(selectedKOmetadata2) == colnames(cts_select2))


##create DEseq dataset
dds3<- DESeqDataSetFromMatrix(countData = cts_select2,
                              colData = selectedKOmetadata2,
                              design = ~ condition)
dds3

##pre-filteringreads <5
keep <- rowSums(counts(dds3))>= 10
dds3<- dds3[keep,]

#specify level of comparison
dds3$condition <- relevel(dds3$condition, ref = "FGSC4200")


#Now, run DEseq and plot dispesion estimates to see the quality of the fit. 

```{r coldata}
dds3 <- DESeq(dds3)
plotDispEsts(dds3)
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/SelectedKO20")
ggsave(filename = "DeSeq_SelectKO2_DispersionEstimate.pdf", plot = dds3, dpi=4200)
```
Export genes with statistically significant differential expression when compared against the wild type.

```{r coldata}
library(DESeq2)
alpha = 0.05

res_crf4 <- results(dds3, alpha=alpha, contrast=c("condition", "crf4" , "FGSC4200"))
res_S238 <- results(dds3, alpha=alpha, contrast=c("condition", "S238" , "FGSC4200"))
res_FGSC12770 <- results(dds3, alpha=alpha, contrast=c("condition", "FGSC12770" , "FGSC4200"))
res_FGSC14852 <- results(dds3, alpha=alpha, contrast=c("condition", "FGSC14852" , "FGSC4200"))


write.csv(as.data.frame(res_crf4), file = "res_crf4.csv")
write.csv(as.data.frame(res_S238), file = "res_S238.csv")
write.csv(as.data.frame(res_FGSC12770), file = "res_FGSC12770.csv")
write.csv(as.data.frame(res_FGSC14852), file = "res_FGSC14852.csv")



list_of_files <- list.files(path = "/Users/ry00555/Desktop/SelectedKO2",
                            recursive = TRUE,
                            pattern = ".csv$")

induced_genes <- readr::read_csv(list_of_files, id = "file_name")


##convert to wide format
l2fc_selectedKO2_data <- data.frame(pivot_wider(data = induced_genes, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_selectedKO2_data) <- l2fc_selectedKO2_data[,1]
l2fc_selectedKO2_data <- l2fc_selectedKO2_data[,-1]
SelectedKO2_PRC2Targetgenesl2fc <-subset(l2fc_selectedKO2_data, rownames(l2fc_selectedKO2_data)%in%K27SilentAcounts[,1])
SelectedKO2_PRC2Targetgenesl2fc <- SelectedKO2_PRC2Targetgenesl2fc[,-(5:8)]

#heat map with log 2 fold change 
selected2 <- read.table("/Users/ry00555/Desktop/SelectedKO2.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
selected_cols <- allDataTPM[,selected2$V1]
GenesWithChanges <- subset(selected_cols, (rowSums(selected_cols) > 0))

heatmap_selected2KO_GenesWChanges<- pheatmap(SelectedKO2_PRC2Targetgenesl2fc, color = colorRampPalette(rev(brewer.pal(n = 10, name="PuOr")))(100), breaks = breaks, cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)
breaks <- c(seq(-20, 20, by=0.4))

heatmap_selected2KO_GenesWChanges_plot<-heatmap_selected2KO_GenesWChanges[[4]]
ggsave(filename = "heatmap_selected2KO_GenesWChanges.pdf", plot = heatmap_selected2KO_GenesWChanges_plot, dpi=4200)

Significanceheatmap_selected2KO_GenesWChanges<- pheatmap(PRC2Targetgenesl2fc, color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100), cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)

ggsave(filename = "Heatmap_Selected.pdf", plot = heatmap_plot3, dpi=4200)


Subset to only export significant differentially expressed genes

```{r}
setwd("/Users/ry00555/Desktop/DeSeqSigSelectedKO2")
dds4<- DESeqDataSetFromMatrix(countData = cts_select2,
                              colData = selectedKOmetadata2,
                              design = ~ condition)
dds4

##pre-filteringreads <5
keep <- rowSums(counts(dds4))>= 10 
dds4<- dds4[keep,]

#specify level of comparison
dds4$condition <- relevel(dds4$condition, ref = "FGSC4200")


#Now, run DEseq and plot dispesion estimates to see the quality of the fit. 

```{r coldata}
dds4 <- DESeq(dds4)
plotDispEsts(dds4)

ggsave(filename = "DeSeq_SelectKO2_DE_Sig.pdf", plot = dds4, dpi=4200)

alpha = 0.05
padj.cutoff <- 0.05
res_sig_crf4 <- subset(res_crf4, padj < 0.05)
res_sig_S238 <- subset(res_S238, padj < 0.05)
res_sig_FGSC12770 <- subset(res_FGSC12770, padj < 0.05)
res_sig_FGSC14852 <- subset(res_FGSC14852, padj < 0.05)


write.csv(as.data.frame(res_sig_crf4), file = "res_sig_crf4.csv")
write.csv(as.data.frame(res_sig_S238), file = "res_sig_S238.csv")
write.csv(as.data.frame(res_sig_FGSC12770), file = "res_sig_FGSC12770.csv")
write.csv(as.data.frame(res_sig_FGSC14852), file = "res_sig_FGSC14852.csv")
list_of_sigfiles <- list.files(path = "/Users/ry00555/Desktop/iswset7",
                               recursive = TRUE,
                               pattern = ".csv$")
setwd("/Users/ry00555/Desktop/iswset7")
induced_genes_sig2 <- readr::read_csv(list_of_sigfiles, id = "file_name")
l2fc_selectedKO2_sigdata <- data.frame(pivot_wider(data = induced_genes_sig2, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_selectedKO2_sigdata) <- l2fc_selectedKO2_sigdata[,1]
l2fc_selectedKO2_sigdata <- l2fc_selectedKO2_sigdata[,-1]
l2fc_selectedKO2_sigdata[is.na(l2fc_selectedKO2_sigdata)]<-0
GenesWithChanges2 <- subset(l2fc_selectedKO2_sigdata, (rowSums(l2fc_selectedKO2_sigdata) > 0))


l2fc_selectedKO2_sigdata <- data.frame(pivot_wider(data = induced_genes_sig2, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_selectedKO2_sigdata) <- l2fc_selectedKO2_sigdata[,1]
l2fc_selectedKO2_sigdata <- l2fc_selectedKO2_sigdata[,-1]
SelectedKO2_PRC2Targetgenes_SIGl2fc <-subset(GenesWithChanges2, rownames(l2fc_selectedKO2_sigdata)%in%K27SilentAcounts[,1])
SelectedKO2_PRC2Targetgenes_SIGl2fc <- SelectedKO2_PRC2Targetgenes_SIGl2fc[,-(5:8)]
SelectedKO2_PRC2Targetgenes_SIGl2fc <- subset(SelectedKO2_PRC2Targetgenes_SIGl2fc, (rowSums(SelectedKO2_PRC2Targetgenes_SIGl2fc) > 0))
library(pheatmap)

#heat map with log 2 fold change 
breaks <- c(seq(0, 10, by=0.1))
SigGenesWChanges <- pheatmap(SelectedKO2_PRC2Targetgenes_SIGl2fc, color = colorRampPalette((brewer.pal(n = 7, name="YlOrRd")))(200), breaks=breaks, cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)

Significanceheatmap_selected2KO_GenesWChanges_plot<-Significanceheatmap_selected2KO_GenesWChanges[[4]]
ggsave(filename = "Significanceheatmap_selected2KO_GenesWChanges.pdf", plot = heatmap_selected2KO_GenesWChanges_plot, dpi=4200)


````````````````````````
11/1/22
Make a MA plot 
lbrary(plotMA)
plotMA(res_sig_crf4, ylim=c(-2,2))
`````````````````````
samplenumbers <-heatmap$tree_col$order[1:20]
samplenames <- heatmap$tree_col$labels[samplenumbers]
write.table(samplenames, file="./Top20.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
SelectedSamplesTop20 <- PRC2TargetTPM[,samplenames]
GenesWithChanges_selectedSamples <- subset(SelectedSamplesTop20, (rowSums(SelectedSamplesTop20) > 0))

````````````````````
setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Rerun11222")
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

Allcountsmatrix <- data.matrix(Allcounts)
Allcounts_countsOnly <- data.matrix(Allcountsmatrix[ ,6:ncol(Allcountsmatrix)])
Allcounts_TPM <- calculateTPM(Allcounts_countsOnly, lengths = Allcountsmatrix[,5])
allDataTPM  <- data.matrix(Allcounts_TPM)
colnames(Allcounts_countsOnly) <- SRRtoGENE[,1]
colnames(allDataTPM) <- SRRtoGENE[,1]


Top20gene <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Top20KOwPRC2GeneChanges.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
Top20geneMeta <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Top20KOwPRC2GeneChangesMetadata.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")

Top20gene_cols <- allDataTPM[,Top20gene$V1]
Top20gene_counts <- Allcounts_countsOnly[,Top20gene$V1]
Top20genechanges <- subset(Top20gene_counts, (rowSums(Top20gene_counts) > 0))


Top20gene_countsmatrix <- as.matrix(Top20genechanges)
Top20gene_numeric <- matrix(as.numeric(Top20genechanges), ncol = 23)
dimnames(Top20gene_numeric) <- list(rownames(Top20genechanges), colnames(Top20genechanges))
Top20_shortened <- as.data.frame(as.matrix(cbind(Top20gene_numeric[,1:3], Top20gene_numeric[,4:5], Top20gene_numeric[,6:8], Top20gene_numeric[,9], Top20gene_numeric[,10:12], Top20gene_numeric[,13:15], Top20gene_numeric[,16:18], Top20gene_numeric[,19:21], Top20gene_numeric[,22], Top20gene_numeric[,23])))

rownames(Top20geneMeta) <- make.names(Top20geneMeta[,1], unique = TRUE)
Top20gene_countsmatrix <- Top20_shortened
colnames(Top20gene_countsmatrix) <- rownames(Top20geneMeta)


Top20geneMeta$condition <- factor(Top20geneMeta$condition)
Top20geneMeta$type <- factor(Top20geneMeta$type)
all(rownames(Top20geneMeta) %in% colnames(Top20gene_countsmatrix))
all(rownames(Top20geneMeta) == colnames(Top20gene_countsmatrix))


dds_top20<- DESeqDataSetFromMatrix(countData = Top20gene_countsmatrix,
                              colData = Top20geneMeta,
                              design = ~ condition)
dds_top20

keep <- rowSums(counts(dds_top20))>= 10
dds_top20<- dds_top20[keep,]

dds_top20$condition <- relevel(dds_top20$condition, ref = "wildtype")

dds_top20 <- DESeq(dds_top20)
plotDispEsts(dds_top20)
ggsave(filename = "Top20_DispEst.pdf", plot = dds_top20, dpi=4200)

alpha = 0.05

res_isw <- results(dds_top20, alpha=alpha, contrast=c("condition", "isw", "wildtype"))
res_ncu06787 <- results(dds_top20, alpha=alpha, contrast=c("condition", "ncu06787", "wildtype"))
res_rtt109 <- results(dds_top20, alpha=alpha, contrast=c("condition", "rtt109", "wildtype"))
res_cac3 <- results(dds_top20, alpha=alpha, contrast=c("condition", "cac3", "wildtype"))
res_H3K4R <- results(dds_top20, alpha=alpha, contrast=c("condition", "H3K4R", "wildtype"))
res_ncu06788 <- results(dds_top20, alpha=alpha, contrast=c("condition", "ncu06788", "wildtype"))
res_ncu041017 <- results(dds_top20, alpha=alpha, contrast=c("condition", "ncu041017", "wildtype"))
res_ncu00548 <- results(dds_top20, alpha=alpha, contrast=c("condition", "ncu00548", "wildtype"))

res_sig_isw <- subset(res_isw, padj < 0.05)
res_sig_ncu06787 <- subset(res_ncu06787, padj < 0.05)
res_sig_rtt109 <- subset(res_rtt109, padj < 0.05)
res_sig_cac3 <- subset(res_cac3, padj < 0.05)
res_sig_H3K4R <- subset(res_H3K4R, padj < 0.05)
res_sig_ncu06788 <- subset(res_ncu06788, padj < 0.05)
res_sig_ncu041017 <- subset(res_ncu041017, padj < 0.05)
res_sig_ncu00548 <- subset(res_ncu00548, padj < 0.05)

write.csv(as.data.frame(res_sig_isw), file = "res_sig_isw.csv")
write.csv(as.data.frame(res_sig_ncu06787), file = "res_sig_ncu06787.csv")
write.csv(as.data.frame(res_sig_rtt109), file = "res_sig_rtt109.csv")
write.csv(as.data.frame(res_sig_cac3), file = "res_sig_cac3.csv")
write.csv(as.data.frame(res_sig_H3K4R), file = "res_sig_H3K4R.csv")
write.csv(as.data.frame(res_sig_ncu06788), file = "res_sig_ncu06788.csv")
write.csv(as.data.frame(res_sig_ncu041017), file = "res_sig_ncu041017.csv")
write.csv(as.data.frame(res_sig_ncu00548), file = "res_sig_ncu00548.csv")

list_of_Top20sigfiles <- list.files(path = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Rerun11222",
                               recursive = TRUE,
                               pattern = ".csv$")

induced_genes_Top20sig <- readr::read_csv(list_of_Top20sigfiles, id = "file_name")
````````````````````````````````
TOP 20 KO STRAINS only show logfold2change in 35 of 679 PRC2 Silenced genes 
l2fc_Top20 <- data.frame(pivot_wider(data = induced_genes_Top20sig, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_Top20) <- l2fc_Top20[,1]
l2fc_Top20 <- l2fc_Top20[,-1]
l2fc_Top20_PRC2Targetonly <-subset(l2fc_Top20, rownames(l2fc_Top20)%in%K27SilentAcounts[,1])
l2fc_Top20_PRC2Targetonly <- subset(l2fc_Top20_PRC2Targetonly, (rowSums(l2fc_Top20_PRC2Targetonly) > 0))

library(WriteXLS)
WriteXLS(l2fc_Top20_PRC2Targetonly, ExcelFileName="/Users/ry00555/Desktop/l2fc_Top20_PRC2Targetonly.xlsx", AdjWidth = T, BoldHeaderRow = T)
l2fc_Top20_PRC2Targetonly.xlsx <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/l2fc_Top20_PRC2Targetonly.xlsx")

tableNA <- write.table(l2fc_Top20_PRC2Targetonly, file="/Users/ry00555/Desktop/l2fc_Top20_PRC2TargetonlyKEEPNA.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
keepNA <- read.table("/Users/ry00555/Desktop/l2fc_Top20_PRC2TargetonlyKEEPNA.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")


breaks <- c(seq(-5, 15, by=0.2))
Top20_PRC2Target_Heatmap<- pheatmap(l2fc_Top20_PRC2Targetonly, color = colorRampPalette(rev(brewer.pal(n = 10, name="RdYlBu")))(200), cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean", legend=T, show_rownames=TRUE, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 5)
Top20_PRC2Target_Heatmaps_plot<-Top20_PRC2Target_Heatmap[[4]]
ggsave(filename = "Top20_PRC2Target_Heatmap.pdf", plot = Top20_PRC2Target_Heatmap, dpi=4200)
````````````````````````````````
#changes the NA to 0 so you can plot all 
#l2fc_Top20_PRC2Targetonly[is.na(l2fc_Top20_PRC2Targetonly)]<-1
`````````````````
PLot log2Fold Changes as a function of expression with significant and non significant read
isw_sig <- read.table(file = "/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/Rerun11222/res_sig_isw.csv")
isw_sig_target_only <- subset(isw_sig, rownames(isw_sig)%in%K27SilentAcounts[,1])
isw_sig_target_only <- subset(isw_sig_target_only, (rowSums(l2fc_Top20_PRC2Targetonly.xlsx) > 0))
alpha = 0.05
plotMA(isw_sig_target_only, alpha = alpha, colSig = "red",colNonSig = "black")
abline(h=c(-2,2))
dev.copy2pdf(file = "res_isw_l2fc_Fexpression.plot")


```

#Try to do a volcano plot

setwd("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/Top20KOPRC2GeneChanges/VolcanoPlotTop20")
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
library(ggrepel)
Allcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JgiAllSampleCounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
K27SilentAcounts <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/JGI_ListofK27SilentGenesNCUs.txt", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
SRRtoGENE <- read.table("/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/JGI/RNASeqAllJGI_KOonly.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")

Allcountsmatrix <- data.matrix(Allcounts)
Allcounts_countsOnly <- data.matrix(Allcountsmatrix[ ,6:ncol(Allcountsmatrix)])
Allcounts_TPM <- calculateTPM(Allcounts_countsOnly, lengths = Allcountsmatrix[,5])
allDataTPM  <- data.matrix(Allcounts_TPM)
colnames(Allcounts_countsOnly) <- SRRtoGENE[,1]
colnames(allDataTPM) <- SRRtoGENE[,1]
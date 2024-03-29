#########R code chunks###################
##R read in bam
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(edgeR)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(stringr)
library(ggbeeswarm)
library(spp)
library(ShortRead)
library(ChIPseeker)
library(forcats)
library(dplyr)

#List of ChIP-Seq Bam files for importing later
#Note, ChIP-Seq libraries were mapped on a merged genome consiting of A.Nidulans and the co-incubated Strapas

pathbase  <- file.path( "/home/sm934/analysis/jule_chipseq" )
pathbam   <- file.path(pathbase,"data")
fls       <- list.files(pathbam,pattern="mapped.bam$",full.names=TRUE)
fls.short <- list.files(pathbam,pattern="mapped.bam$",full.names=FALSE)

##reading in genome, annotation.
genome           <- readDNAStringSet(file.path(pathbase,"/Anidulans_Siranensis_Srapamucinicus_merged_genomes.fa"))
gffgenes         <- import.gff3(file.path(pathbase,"/A_nidulans_FGSC_A4_version_s10-m03-r28_features.gff"))
chr_sizes        <- width(genome)
names(chr_sizes) <- sapply(strsplit(names(genome)," "),function(x) x[[1]])
mRNA <- gffgenes[gffgenes@elementMetadata$type=="mRNA",]
tes  <- mRNA[grep("transposon",unlist(mRNA$Note)),]
downstream <- 500
upstream <- 1000
promotors5do10up <- mRNA
promotors5do10up $type="promotor"
for (i in 1:length(promotors5do10up)){
  if(as.vector(strand(promotors5do10up[i])=="+")) {
    TSS <- start(promotors5do10up[i])
    start(promotors5do10up[i]) <- TSS - upstream
    end(promotors5do10up[i]) <- TSS + downstream
  } else {
    TSS <- end(promotors5do10up[i])
    start(promotors5do10up[i]) <- TSS - downstream
    end(promotors5do10up[i]) <- TSS + upstream
  }
}
#Note, there is a promotor function in Biostrings or so, have to redo next time
load("Robjects_nidulans.rdata")
#save(chr_sizes,genome,fls,gffgenes,febitA,febitF,chip_merged,chiplogFCsall,cnts,probenames,counts_edger,countsnall,groupdf,groupdf2,gwcov, cov.5e4,cov.5e5,cov.2e3,cov.5e4.norm,cov.2e3.norm,cov.5e5.norm,file="Robjects_nidulans.rdata")

# reading in stats table
stats_counts <- read.table(file.path(pathbase,"data/Bowtie2_mapping_rates_chipseq_anidulans.csv"),header=T)
stats_mito <- 100*stats_counts[grep("mito",stats_counts[,1]),-1] / colSums(stats_counts[,-1])
#   GcnE.FLAG_noSrapa_noSiran_H3Cterm_Rep1 GcnE.FLAG_noSrapa_noSiran_H3Cterm_Rep2
#                               0.0179859                             0.01362185
#  GcnE.FLAG_noSrapa_noSiran_H3Cterm_Rep3 GcnE.FLAG_noSrapa_noSiran_H3K14ac_Rep1
#                              0.01553365                             0.03460474
#  GcnE.FLAG_noSrapa_noSiran_H3K14ac_Rep2 GcnE.FLAG_noSrapa_noSiran_H3K14ac_Rep3
#                              0.02910052                             0.05029246
#  GcnE.FLAG_noSrapa_noSiran_H3K9ac_Rep1 GcnE.FLAG_noSrapa_noSiran_H3K9ac_Rep2
#                             0.03786612                            0.03447385
#  GcnE.FLAG_noSrapa_noSiran_H3K9ac_Rep3 GcnE.FLAG_Srapa_noSiran_FLAGM2_Rep2
#                             0.01741689                           0.2615729
#  GcnE.FLAG_Srapa_noSiran_H3Cterm_Rep1 GcnE.FLAG_Srapa_noSiran_H3Cterm_Rep2
#                            0.01227476                            0.0122743
#  GcnE.FLAG_Srapa_noSiran_H3Cterm_Rep3 GcnE.FLAG_Srapa_noSiran_H3K14ac_Rep1
#                           0.009824139                           0.03225528
#  GcnE.FLAG_Srapa_noSiran_H3K14ac_Rep2 GcnE.FLAG_Srapa_noSiran_H3K14ac_Rep3
#                            0.01384039                           0.02542382
#  GcnE.FLAG_Srapa_noSiran_H3K9ac_Rep1 GcnE.FLAG_Srapa_noSiran_H3K9ac_Rep2
#                            0.0361413                          0.02798489
#  GcnE.FLAG_Srapa_noSiran_H3K9ac_Rep3 WT_Srapa_noSiran_FLAGM2_Rep1
#                            0.0306326                    0.2016542
#  WT_Srapa_noSiran_FLAGM2_Rep1.1
#                       0.2586124
# from ~0.02% (H3) to ~0.2% (FlagM2)
# ~ 10 fold enrichment control(flag) vs H3 etc

#some QC with SPP : cross correlation plots
#http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html
# require(devtools)
# devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)
#this needs to be restricted to only nidulans genome for proper estimation
peak_sep <- rep(NA, length(BamFileList(fls)))
for (i in 1:length(BamFileList(fls))) {
  print(BamFileList(fls)[[i]]$path)
  name <- rev(strsplit(BamFileList(fls)[[i]]$path,"/")[[1]])[1]
  print(name)
  alnspp <- read.bam.tags(BamFileList(fls)[[i]]$path, read.tag.names = F, fix.chromosome.names = F)
  binding.characteristics <- get.binding.characteristics(alnspp,srange=c(50,500),bin=5)
  # print out binding peak separation distance
  print(paste("binding peak separation distance =",binding.characteristics$peak$x))
  peak_sep[i] <- binding.characteristics$peak$x
  # plot cross-correlation profile
  png(file=paste0(name,".png"),width=700,height=500)
  par( mgp = c(2,0.65,0), cex = 0.8)
  plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation",sub=name)
  abline(v=binding.characteristics$peak$x,lty=2,col=2)
  dev.off()
}


#####R analysis for DiffBinding#########
###reading in BAM files for counting####

#testing parameters to read in bam files, 
paramtest               <- ScanBamParam(flag=scanBamFlag(isNotPrimaryRead=FALSE, isProperPair=T,isDuplicate=F,isNotPassingQualityControls=F,isMinusStrand=NA,isFirstMateRead=NA),tag=c("NH", "AS", "NM"),what=c("mapq", "flag"))
aln                     <- readGAlignmentsFromBam(BamFileList(fls)[[1]], use.names=T, paramtest=paramtest)
alnspp                  <- read.bam.tags(BamFileList(fls)[[1]]$path, read.tag.names = F, fix.chromosome.names = F)
binding.characteristics <- get.binding.characteristics(alnspp,srange=c(50,500),bin=5)
##select testing intervals (plus a name to specify file names), run code below 
##needs $Name slot 

param <- ScanBamParam(flag=scanBamFlag(isDuplicate=F,isNotPassingQualityControls=F,isMinusStrand=NA,isFirstMateRead=NA),tag=c("NH", "AS", "NM"),what=c("mapq", "flag"))

#We are testing wheater counting is best done on Genobodies or promoters
#the same anlysis can be performed for both and can be changed by setting "label" to Genes or Promotors:
intervals <- mRNA; label <- "Genes" ##set globaly

gnCntgenes <- summarizeOverlaps(intervals, BamFileList(fls), mode="Union",ignore.strand=TRUE, singleEnd=TRUE,param=param)
intervals  <- promotors5do10up
gnCntproms <- summarizeOverlaps(intervals, BamFileList(fls), mode="Union",ignore.strand=TRUE, singleEnd=TRUE,param=param)
save(gnCntgenes,gnCntproms,file="gene_prom_counts.rdata")

load("gene_prom_counts.rdata")

if(label=="Genes") {gnCnt <- gnCntgenes} else {gnCnt <- gnCntproms}

#renaming and extracting counts
namestmp <- as.character(1:length(intervals))
namestmp[!is.na(intervals$Name)] <- intervals$Name[!is.na(intervals$Name)]
rownames(gnCnt) <- namestmp
cnts <- assay(gnCnt) #gnCnt@assays$data@listData$counts
cnts <- cnts[,-1] #excluding first sample since it is redundant# 
dim(cnts) # [1] 10771    43 

#setting up metadata and design structure for experiment
groupnames                          <- strsplit(colnames(cnts),"_")
groupdf                             <- NULL
groupdf$strain                      <- as.factor(sapply(groupnames,function(x) x[[1]]))
groupdf$S_Rapa_coincubation         <- as.factor(sapply(groupnames,function(x) x[[2]]))
levels(groupdf$S_Rapa_coincubation) <- c("-S.rapa","+S.rapa")
groupdf$S_Iran_coincubation         <- as.factor(sapply(groupnames,function(x) x[[3]]))
levels(groupdf$S_Iran_coincubation) <- c("-S.iran","+S.iran","+S.iran-mutant")
groupdf$Antibody                    <- as.factor(sapply(groupnames,function(x) x[[4]]))
groupdf$Replicate                   <- as.factor(sapply(groupnames,function(x) x[[5]]))
groupdf$ID                          <- sapply(groupnames,function(x) x[[6]])
groupdf$condition                   <- as.factor(paste(groupdf$S_Iran_coincubation,groupdf$S_Rapa_coincubation))
levels(groupdf$condition)           <- c("GcnE-3xFLAG", "GcnE-3xFLAG + S.rapamycinicus","GcnE-3xFLAG + S.iranensis", "GcnE-3xFLAG + S. iranensis mutant")
groupdf           <- as.data.frame(groupdf) #43
groupdf$libnames  <- colnames(cnts)
# groupdf$labels <- groupdf$condition
# levels(groupdf$condition)
#taking out FLAGM2 since it didn't work well
subsetFlag       <- groupdf$Antibody != "FLAGM2" #43
groupdfNoFlag    <- subset(groupdf, subsetFlag)
dim(groupdfNoFlag) # [1] 36  7
#       strain S_Rapa_coincubation S_Iran_coincubation Antibody Replicate    ID                         condition     gene
# 1  GcnE-FLAG             -S.rapa             -S.iran  H3Cterm      Rep1 25985                       GcnE-3xFLAG 487.9121


###DCB analysis using edgeR######
group <- as.factor(sapply(strsplit(colnames(cnts),"_Rep"),function(x) x[[1]]))
counts_edger <- DGEList(counts=cnts,group=group)
# cntsnc  =  cpm(cnts_edger, normalized.lib.sizes = TRUE)
counts_edger <- calcNormFactors(counts_edger)
counts_edger <- estimateCommonDisp(counts_edger)
counts_edger <- estimateTagwiseDisp(counts_edger)

# NORMALIZING
# note, in contrast to the normalizing stratgy for the coverage, Strepto reads (i.e. undesired mapped reads) don't need to be excluded explicitly since the DGEList only contains genes from A. Nidulans which indirectly filters it. 
# NormFactor is therefore calculated only taking A.Nidulans mapped reads into account
lib.sizes <- counts_edger$sample$lib.size * calcNormFactors(counts_edger)$samples$norm.factors
countsnall <- t(t(cnts * mean(lib.sizes))  / lib.sizes)

countsNoFlag  <- countsnall[,subsetFlag] #36
# filtering out libraries with S. iranensis (not used for paper this paper)
subsetSiran <- grep("_Siran",colnames(countsNoFlag)) #  [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
namessubsetSiran <- colnames(countsNoFlag)[-subsetSiran] # 18

##normalzing reads (TMM from edgeR)
#dividing each column of the count table by the correspond-ing size factor yields normalized count values, which can be scaled to give a counts per million interpretation (see also edgeR’s cpm function)

compare_strings_mat <- cbind(  c( "_noSrapa_noSiran_","_Srapa_noSiran_"),
                               c( "_noSrapa_Siran_","_noSrapa_Siran-mutant_"),
                               c( "_noSrapa_noSiran_","_noSrapa_Siran_"),
                               c( "_noSrapa_noSiran_","_noSrapa_Siran-mutant_")
)
myfdr <- 0.01
mylfc <- 0
##comparing A.Nidulans with and without S.rapa for K9, K14 and H3 (as control)
chip_comparisons <- list()
#first element in vector is without S.rapa second with S.rapa 
compare_strings <- c("_noSrapa_noSiran_","_Srapa_noSiran_")
for (AB in c("H3K9ac","H3K14ac","H3Cterm")) {
  countsab <- cnts[,grep(AB,colnames(cnts))]
  #	for (i in 3:ncol(compare_strings_mat)) {
  #		compare_strings <- compare_strings_mat[,i]
  #counts <- countsab[,grep("noSiran",colnames(countsab))]
  counts <- countsab[,c(grep(compare_strings[1],colnames(countsab)),grep(compare_strings[2],colnames(countsab)))]
  #group <- as.factor(rep(paste0(c("noSrapa_noSiran_","Srapa_noSiran_"),AB),each=3))
  group <- as.factor(sapply(strsplit(colnames(counts),"_Rep"),function(x) x[[1]]))
  counts_edger <- DGEList(counts=counts,group=group)
  # cntsnc  =  cpm(cnts_edger, normalized.lib.sizes = TRUE)
  counts_edger <- calcNormFactors(counts_edger)
  counts_edger <- estimateCommonDisp(counts_edger)
  counts_edger <- estimateTagwiseDisp(counts_edger)
  et <- exactTest(counts_edger)
  et$table$decideTestsDGE <- as.vector(decideTestsDGE(et, p.value=myfdr,adjust.method="fdr",lfc=mylfc)) #tags (up/no/down)
  mygenes <- factor(et$table$decideTestsDGE,levels=c(-1,0,1)) %>%
    fct_recode("blue"="-1","grey"="0","red"="1")
  chip_comparisons[[AB]] <- et
  #MAplot
  png(file=paste0("MAplot_norm",AB,".png"),width=1000,height=1000)
  par(mfrow=c(1,1),mar=c(7.1,7.1,4.1,2.1))
  levels(counts_edger$samples$group) <- c("S.rapa(+) "," S.rapa(-)")
  # [1] "GcnE-FLAG_noSrapa_noSiran_H3Cterm" "GcnE-FLAG_Srapa_noSiran_H3Cterm"  
  message(AB," logFC: ",mylfc," FDR: ", myfdr)
  print(table(mygenes))
  # max y range was found to be logFC=2.5: using for consistent axis (ylim)
  plotSmear(counts_edger,lowess=T,cex.lab=2, cex.axis=2, main=AB,cex.main=3,col=as.character(mygenes),ylim=c(-2.5,2.5))
  abline(h = 0, col = "red")
  abline(h = c(-1, 1), col = "blue")
  dev.off()
  #		lib.sizes <- counts_edger$sample$lib.size * calcNormFactors(counts_edger)$samples$norm.factors
  #		countsn <- t(t(counts * mean(lib.sizes))  / lib.sizes)
  #		tmp <- data.frame(topTags(et,sort.by="none",n=length(intervals)),countsn,as.data.frame(intervals))[order(et$table$PValue),]
  #		write.csv2(tmp,file=paste0("diffbinding_edgeR_",AB,paste(compare_strings,sep = "",collapse="vs"),label,".csv"))
  #browser()
  #		pdf(file=paste0("diffbinding_edgeR_",AB,paste(compare_strings,sep = "",collapse="vs"),label,"_pvalues.pdf"))
  #		hist(et$table$PValue,breaks=300,main=paste0("Differential-Binding\n",paste(compare_strings,sep = "",collapse=" vs\n")),xlab=paste0("P-Value distribution of",AB," for ", label))
  #		dev.off()
  #	}
}
# H3K9ac logFC: 0 FDR: 0.01
# mygenes
# blue grey  red
#  297 9881  593
# H3K14ac logFC: 0 FDR: 0.01
# mygenes
#  blue  grey   red
#   154 10513   104
# H3Cterm logFC: 0 FDR: 0.01
# mygenes
#  blue  grey   red
#    37 10732     2
# comparing with logFC=1 and FDR=0.1
# H3K9ac logFC: 1 FDR: 0.1
# mygenes
#  blue  grey   red
#    50 10568   153
# H3K14ac logFC: 1 FDR: 0.1
# mygenes
#  blue  grey   red
#     0 10771     0
# H3Cterm logFC: 1 FDR: 0.1
# mygenes
#  blue  grey   red
#     6 10765     0

countsn_edger <- DGEList(counts=countsn,group=groups)
png(file=paste0("comparison_readcounts_boxplot_noramlization",mapping,".png"))
par(mfcol=c(1:2))
boxplot(log(counts+1),main="raw read log2 counts")
boxplot(log(countsn+1),main="TMM normalized read log2 counts")
dev.off()
#cnts
png(file=paste0("mds_rnaseqsamples.png",".png"),width=500,height=500)
plotMDS(counts, labels = groups,)
dev.off()
png(file=paste0("mds_rnaseqsamples_norm.png",".png"),width=500,height=500)
plotMDS(countsn, labels = groups,)
dev.off()
png(file=paste0("MAplot_norm.png",".png"),width=3000,height=3000)
par(mfcol=c(5,5),cex.lab=2)
for(i in 1:5) {
  for(j in 1:5) {
    if (i==j) { hist(log2(countsn[,i]),breaks=200,sub=colnames(counts)[i]) } else {
      plotSmear(countsn_edger[,c(i,j)],lowess=TRUE,sub=paste(colnames(counts)[i],"vs",colnames(counts)[j]))
      abline(h = 0, col = "red")
      abline(h = c(-2, 2), col = "blue")
    }
  }
}
dev.off()

###Microarray data from Nuetzmann 2011
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25266
febitA     <- read.csv(file.path(pathbase,"/febit/FEAS020/results/comparisons_comparisonA.csv"))
febitF     <- read.csv(file.path(pathbase,"/febit/FEAS020/results/comparisons_comparisonF.csv"))
probenames <- read.csv(file.path(pathbase,"/febit/FEAS020/febit-nidulans-2009_bg_GENID.txt",sep="\t"))
febitF$gene.name <- febitA$gene.name <- sapply(strsplit(as.character(probenames$X),"\\."),function(x) x[[1]])

geneidx <- which(febitF$gene.name %in% genelist)
write.csv(febitF[which(febitF$gene.name %in% genelist),],file="febit_comparisonF_testlist.csv")

plot(febitF$logqmedian,febitA$logqmedian)
text(febitF$logqmedian[geneidx],febitA$logqmedian[geneidx],genelist,col="red")

tmp2          <- febitA[which(febitA$gene.name %in% genelist),]
tmp2          <- febitF[which(febitF$gene.name %in% genelist),]
tmp           <- et$table[genelist,]
tmp$gene.name <- rownames(tmp)
tmp3          <- merge(tmp,tmp2, by="gene.name")
plot(tmp3$logFC,-tmp3$logqmedian,col="grey")
text(tmp3$logFC,-tmp3$logqmedian,tmp3$gene.name)
cor(tmp3$logFC,-tmp3$logqmedian)

##correlation between chip-seq and rna-seq
plot_list <- list()
for (AB in c("H3K9ac","H3K14ac","H3Cterm")) {
  chiptmp           <- chip_comparisons[[AB]]$table
  chiptmp$gene.name <- rownames(chiptmp)
  chip_merged       <- merge(chiptmp,febitF, by="gene.name")
  #febitF$limma_adjp
  tmp            <- chip_merged[,c("logFC","logqmedian")]
  tmp$logqmedian <- -tmp$logqmedian
  rownames(tmp)  <- chip_merged$gene.name
  pdf(file=paste("febit_vs_chipseq_foldchange_",AB,"_",label,".pdf",sep=""),width=12, height=12)
  plot(tmp,col=(rownames(tmp)%in% genelist)+1)
  text(tmp[genelist,],rownames(tmp[genelist,]),pos=1)
  text(tmp[which(abs(tmp$logqmedian) > 6),],rownames(tmp[which(abs(tmp$logqmedian) > 6),]),col="blue",pos=1)
  dev.off()
  message(AB)
  #     print(cor((tmp[chip_merged$PValue<0.001,]))) #0.39
  message("all genes")
  print(cor(tmp[chip_merged$limma_adjp<00.1,],method="spearman")) #0.45
  message("DCB genes")
  print(cor(tmp),method="spearman") #0.45
  #0.56 for febitA instead febitF
}
# H3K9ac
# all genes
#                logFC logqmedian
# logFC      1.0000000  0.4969893
# logqmedian 0.4969893  1.0000000
# DCB genes
#               logFC logqmedian
# logFC      1.000000   0.239166
# logqmedian 0.239166   1.000000
# H3K14ac
# all genes
#                  logFC  logqmedian
# logFC       1.00000000 -0.04821228
# logqmedian -0.04821228  1.00000000
# DCB genes
#                  logFC  logqmedian
# logFC       1.00000000 -0.04149645
# logqmedian -0.04149645  1.00000000
# H3Cterm
# all genes
#                 logFC logqmedian
# logFC      1.00000000 0.05431411
# logqmedian 0.05431411 1.00000000
# DCB genes
#                  logFC  logqmedian
# logFC       1.00000000 -0.02415712
# logqmedian -0.02415712  1.00000000
# 


chiplogFCsall <- data.frame(H3K9ac=chip_comparisons[["H3K9ac"]]$table$logFC,H3K14ac=chip_comparisons[["H3K14ac"]]$table$logFC,H3Cterm=chip_comparisons[["H3Cterm"]]$table$logFC)
chiplogFCsall$gene.name<-rownames(chip_comparisons[["H3K9ac"]]$table)
chip_merged <- merge(chiplogFCsall,febitF[,c("logqmedian","gene.name")], by="gene.name")
chip_merged$logqmedian <- -chip_merged$logqmedian

## MA plot for microarray to check if it has been normalized. Using median.g1 and median.g2 columns
countsn_edger <- DGEList(counts=chip_merged[,c("median.g1","median.g2")],group=c("median.g1","median.g2"))

png(file=paste0("MAplot_norm_febit",".png"),width=1000,height=1000)
edgeR::plotSmear(countsn_edger,lowess=TRUE)
abline(h = 0, col = "red")
abline(h = c(-2, 2), col = "blue")
dev.off()

chip_merged_wide <- melt(chip_merged,id.vars=c("gene.name","logqmedian"))

genelist2 <- c("AN7909","AN12004","AN7911","AN7912","AN7913","AN7914")
gg <- ggplot(chip_merged_wide[chip_merged_wide$gene.name%in%genelist2,],aes(x=value,y=logqmedian)) +
  geom_point(data=chip_merged_wide, color=mycols["grau"]) +
  geom_point(color=mycols["green"],size=2) +
  geom_text_repel(aes(label = gene.name),box.padding = unit(0.75, "lines"),force=6) +
  labs(x = "log fold change (WT/WT+S.Rapa) ChIP-Seq", y = "log fold change (WT/WT+S.Rapa) Microarray") +
  facet_grid(.~variable)
save_plot(file="test_logFC_microarray_vs_chipseq.pdf",gg,base_width=10,base_height=4)
# geom_hline(yintercept=c(-2,0,2),color="grey")

###
######profile plotsk######
proms_mRNA_3k3k     <- promoters(mRNA, upstream=3000, downstream=3001)
mRNA_tts            <- mRNA
start(mRNA_tts)     <- end(mRNA)
mRNA_tts2           <- resize(mRNA,1,fix="end")
proms_mRNA_3k3k_tts <- promoters(mRNA_tts2, upstream=3000, downstream=3001)

tagMatrixlist_reads <- list()
for (i in fls.short){
  message(i)
  tagMatrixlist_reads[[i]] <- getTagMatrix(granges(readGAlignments(i)),windows=proms_mRNA_3k3k)
}
tagMatrixlist_reads_tts <- list()
for (i in fls.short){
  message(i)
  tagMatrixlist_reads_tts[[i]] <- getTagMatrix(granges(readGAlignments(i)),windows=proms_mRNA_3k3k_tts)
}
##subseting taking out Sir
tagMatrixlist_reads_subset <- tagMatrixlist_reads[namessubsetSiran]
tagMatrixlist_reads_tts_subset <- tagMatrixlist_reads_tts[namessubsetSiran]

ggtss <- plotAvgProf(tagMatrixlist_reads_subset[1:9], xlim=c(-3000, 3000),  xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#ggsave(file="tss_profiles_without_interaction.pdf")
ggtts <- plotAvgProf(tagMatrixlist_reads_tts_subset[1:9], xlim=c(-3000, 3000),  xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#ggsave(file="tts_profiles_without_interaction.pdf")

tagMatrix_H3K14ac = Reduce(function(x,y) x + y, tagMatrixlist_reads[grep("H3K14ac",names(tagMatrixlist_reads))],matrix(0,nrow=10756,ncol=6001))
tagMatrix_H3K9ac = Reduce(function(x,y) x + y, tagMatrixlist_reads[grep("H3K9ac",names(tagMatrixlist_reads))],matrix(0,nrow=10756,ncol=6001))
plotAvgProf(tagMatrix_H3K14ac, xlim=c(-3000, 3000),  xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")> plotAvgProf(tagMatrix_H3K9ac, xlim=c(-3000, 3000),  xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
tagMatrix_H3 = Reduce(function(x,y) x + y, tagMatrixlist_reads_subset[grep("H3Cterm",names(tagMatrixlist_reads))],matrix(0,nrow=10756,ncol=6001))

mycols <- c(blau="#253494",green="#238B45",grau="#969696",schwarz="#000000")
mythemecols <- scale_colour_manual(values=rev(rep(as.character(mycols[1:3]),each=3)),
                                   name="Experimental\nCondition",
                                   breaks=names(tagMatrixlist_reads_tts_subset)[1:9],
                                   labels=rep(c("H3","H3K14ac","H3K9ac"),each=3))

mythemeTTS <- scale_x_continuous(breaks=c(-3000,-2000,-1000,0,1000,2000,3000),labels=c("-3kb","-2kb","-1kb","TTS","1kb","2kb","3kb")) 
mythemeTSS <- scale_x_continuous(breaks=c(-3000,-2000,-1000,0,1000,2000,3000),labels=c("-3kb","-2kb","-1kb","TSS","1kb","2kb","3kb")) 

ggtts <- ggtts + mythemecols + mythemeTTS + theme(legend.position=c(.8,.7)) + ylim(0.00012,0.000265)
ggtss <- ggtss + mythemecols + mythemeTSS + theme(legend.position="none")+ ylim(0.00012,0.000265)
gg <- plot_grid(ggtss, ggtts, labels = c("A", "B"), align = "h")
save_plot(file="profile_plots_tss_tts.pdf",gg,base_width=10,base_height=4)

######coverage#######
bins2e3      <- tileGenome(chr_sizes, tilewidth=2e3, cut.last.tile.in.chrom=T)
bins2e3$ID   <- paste(seqnames(bins2e3),start(bins2e3),sep="_")
bins2e3$type <- "Bin2e3"
bins5e4      <- tileGenome(chr_sizes, tilewidth=5e4, cut.last.tile.in.chrom=T)
bins5e4$ID   <- paste(seqnames(bins5e4),start(bins5e4),sep="_")
bins5e4$type <- "Bin5e4"
bins5e5      <- tileGenome(chr_sizes, tilewidth=5e5, cut.last.tile.in.chrom=T)
bins5e5$ID   <- paste(seqnames(bins5e5),start(bins5e5),sep="_")
bins5e5$type <- "Bin5e5"

#export.gff3(promotors5do10up,con=paste0("promotors5do10up_AN_upstream_",upstream,"_downstream_",downstream,".gff3"))
##Genomewide coverage using bins
cov.5e4 <- gwcov(bins5e4,groupdf)
cov.5e5 <- gwcov(bins5e5,groupdf)
cov.2e3 <- gwcov(bins2e3,groupdf)

##normalizing
￼
#cluster: II:96009-225088
#lib.sizes <- colSums(cov.5e4$wide[,1:nrow(groupdf)])
#only using Nidulans mapped 
# this the key normalizing step which excludes reads mapped to the Strepto!!
lib.sizes <- colSums(cov.5e4$wide[-grep("Streptomyces",cov.5e4$wide$seqnames),1:nrow(groupdf)])

#library size normalizing for coverage to account for differnet sequencing depths 
cov.5e4.norm <- cov.5e4
cov.5e4.norm$wide[,1:nrow(groupdf)] <- t(t(cov.5e4$wide[,1:nrow(groupdf)] * mean(lib.sizes))  / lib.sizes)
cov.5e4.norm$long <- gwcov.long(cov.5e4.norm$wide,groupdf)

cov.2e3.norm <- cov.2e3
cov.2e3.norm$wide[,1:nrow(groupdf)] <- t(t(cov.2e3$wide[,1:nrow(groupdf)] * mean(lib.sizes))  / lib.sizes)
cov.2e3.norm$long <- gwcov.long(cov.2e3.norm$wide,groupdf)

#genome wide coverage plots
ggplot(subset(cov.5e4.norm$long, seqnames %in% c("ChrIII_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(Replicate))) + geom_line() + facet_grid(condition~Antibody)

gg <- ggplot(subset(cov.5e4.norm$long, seqnames %in% c("ChrI_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(Replicate))) + geom_line() + facet_grid(condition~Antibody)
ggsave(gg,file="Genomewide_coverage_5e4_bins_Chrom1_nidulans_normalized.png",width=14)

gg <- ggplot(subset(cov.5e4.norm$long,! seqnames %in% c("mito_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(Replicate))) + geom_line() + facet_grid(condition+Antibody~seqnames)
ggsave(gg,file="Genomewide_coverage_5e4_bins_allChroms_nidulans_normalized.png",width=30,height=30)

gg <- ggplot(subset(cov.5e4.norm$long, seqnames %in% c("ChrII_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(condition))) + geom_smooth(span=0.05) + facet_grid(.~Antibody)
ggsave(gg,file="Genomewide_coverage_smooth_5e4_bins_Chrom2_nidulans_normalized.png",width=30,height=10)

gg <- ggplot(subset(cov.2e3.norm$long,! seqnames %in% c("mito_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(condition))) + geom_smooth(method="loess",span=0.05) + facet_grid(Antibody~.)+ geom_vline(xintercept = c(207766,214402))
ggsave(gg,file="Genomewide_coverage_smooth_2e3_bins_Chrom2_nidulans_normalized.png",width=12,height=4)

gg <- ggplot(subset(cov.2e3.norm$long, seqnames %in% c("ChrII_A_nidulans_FGSC_A4")),aes(y=value,x=start,color=as.factor(condition))) + geom_line() + facet_grid(Antibody~.)+ geom_vline(xintercept = c(207766,214402))

ggsave(gg,file="Genomewide_coverage_line_2e3_bins_Chrom2_nidulans_normalized.png",width=40,height=4)

##Plotting quanititative ChIP signal (to confirm) 
#selected by Juliane based on diff Binding analysis
genelist <- c("AN8772", "AN11891", "AN3675", "AN1812", "AN7944", "AN2161", "AN6411", "AN1899", "AN0992", "AN7327", "AN3125", "AN1160", "AN8489", "AN3124", "AN3587", "AN1330", "AN2505", "AN8102", "AN2903", "AN2919", "AN6642", "AN6132", "AN1008", "AN1006", "AN1007", "AN2944", "AN5134", "AN4376", "AN5696", "AN7463", "AN0232","AN1731", "AN1732", "AN4727", "AN1137", "AN9506", "AN1923", "AN5731", "AN6255", "AN8118", "AN2549", "AN10080", "AN5520", "AN4073", "AN7913", "AN12004", "AN7911", "AN7912", "AN7914", "AN7909", "AN7907", "AN7883", "AN7895", "AN7897", "AN7896")
# genelist <-paste0("AN79",c("09","11","12","13","14")) #ors Genes for Fig2

# genes belong to BGC or and nitron/calcim pathways
genelistall <- list(
  nitrogen_fig5 =  paste0("AN",c("1006","1007","1008","1731","1732","2944","4376","5134","7463")), #other Genes for Fig5
  # biosynthetic
  calcium_signal=  paste0("AN",c("3585","4920","3998","3420","4418","2427","0419","3751","5372","5993","5341","2826","1950","5302","8842","7464","8774")),
  ors_fig2      =  paste0("AN79",c("09","11","12","13","14")), #genes in orsellinic BGC
  atn           =  paste0("AN78",95:97), #atn Genes
  cichorin      =  paste0("AN64",c("27","40","41","47")), #ors
  stc           =  paste0("AN78",c("11","16","21")),
  mdp           =  paste0("AN0147"),
  eas           =  paste0("AN25",48:49)
)

## plotting individual genes of interest:
groupdfTmp <- groupdfNoFlag
for (genename in genelist) {
  groupdfTmp$gene <- countsNoFlag[genename,]
  groupdf2 <- groupdfTmp[-subsetSiran,]
  message(genename)
  #	gg <- ggplot(groupdf, aes(fill=strain,y=gene,x=S_Rapa_coincubation)) + geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",binwidth = 1,dotsize=30) + facet_grid(Antibody~S_Iran_coincubation) + labs(title = genename)
  #	
  #	gg <- ggplot(groupdf, aes(fill=S_Rapa_coincubation,y=gene,x=S_Rapa_coincubation)) + stat_summary(fun.y="mean", geom="bar") +geom_jitter() + facet_grid(Antibody~S_Iran_coincubation) + labs(title = genename,x = "",y="normalized read counts")
  #	
  #	gg <- ggplot(groupdf, aes(fill=condition,y=gene,x=condition)) + stat_summary(fun.y="mean", geom="bar") + facet_grid(Antibody~.) + geom_jitter() +  labs(title = genename,x = "",y="normalized read counts") + scale_fill_manual(values=brewer.pal(4,"RdYlBu"))+  theme_bw()
  #	
  #gg <- ggplot(groupdf, aes(fill=condition,y=gene,x=condition)) + stat_summary(fun.y="mean", geom="bar") + facet_grid(Antibody~.) + geom_quasirandom() +  labs(title = genename,x = "",y="normalized read counts") + scale_fill_manual(values=rev(brewer.pal(4,"RdYlBu"))) + theme_bw() + theme(legend.position="none",legend.direction = "vertical",legend.position="top",axis.text.x = element_text(angle = 45,vjust=0.5)) + scale_x_discrete(labels=gsub("\\+","\n\\+",levels(groupdf$condition)))
  gg <- ggplot(groupdf2, aes(fill=condition,y=gene,x=condition)) +
    stat_summary(fun.y="mean", geom="bar") +
    facet_grid(Antibody~.) +
    geom_jitter(position = position_jitter(width = .03),size=1,alpha=0.7) +
    labs(title = genename,x = "",y="normalized read counts") +
    scale_fill_manual(values=rev(brewer.pal(4,"RdYlBu"))) +
    theme_bw() +
    theme(legend.direction = "vertical",legend.position="top",axis.text.x = element_text(angle = 45,vjust=0.5)) +
    scale_x_discrete(labels=gsub("\\+","\n\\+",levels(groupdf2$condition)))
  #browser()
  #ggsave(gg,file=paste0("read_counts_",genename,"_",label,".pdf"),width=10, height=8)
  ggsave(gg,file=paste0("read_counts_",genename,"_",label,".pdf"),width=2, height=5.5)
}

#create and plot compsite barplot for the 5 ors Genes
# more general now: for all cluster genes (plus nitro and calcium)
for (genelistname in names(genelistall)) {
  #genelistname <- names(genelistall)[1]
  genelist <- genelistall[[genelistname]]
  groupdf2 <- cbind(t(countsNoFlag[genelist,-subsetSiran]),groupdfNoFlag[-subsetSiran,])
  groupdf2$gene <- NULL
  # groupdf2long <- dplyr::rename(groupdf2, orsA = AN7909, orsB = AN7911, orsC = AN7912, orsD = AN7913, orsE = AN7914) %>%
  # groupdf2long <- dplyr::rename(groupdf2, niaD = AN1006, niiA = AN1007, crnA = AN1008, prnD = AN1731, prnB = AN1732, tamA = AN2944, gdhA = AN4376, gltA = AN5134, meaA = AN7463) %>%
  groupdf2long <- groupdf2 %>% 
    gather("Gene","Counts",1:length(genelist)) %>%
    mutate(Antibody = str_replace(Antibody,"Cterm","")) %>%
    mutate(condition = str_replace(condition,"GcnE-3xFLAG","A. nidulans"))
  #adding microarray data
  groupdf2mic <- filter(chip_merged,gene.name %in% genelist) %>%
    # groupdf2mic <- filter(chip_merged,gene.name %in% sample(gene.name,15)) %>%
    select(gene.name,median.g1,median.g2) %>%
    gather("condition","Counts",2:3) %>%
    #   mutate(Counts = log2(Counts)) %>%
    mutate(Replicate = "Rep1") %>%
    #   mutate(strain = "GcnE-FLAG") %>%
    mutate(Antibody = "Microarray") %>%
    mutate(Gene=gene.name) %>%
    #     mutate(Gene=fct_recode(gene.name, orsA = "AN7909", orsB = "AN7911", orsC = "AN7912", orsD = "AN7913", orsE = "AN7914")) %>%
    #   mutate(Gene=fct_recode(gene.name, niaD = "AN1006", niiA = "AN1007", crnA = "AN1008", prnD = "AN1731", prnB = "AN1732", tamA = "AN2944", gdhA = "AN4376", gltA = "AN5134", meaA = "AN7463")) %>%
    mutate(condition=fct_recode(condition,"A. nidulans" = "median.g1","A. nidulans + S.rapamycinicus" = "median.g2")) %>%
    select(-gene.name)
  groupdf2comb <- rbind(groupdf2long[,colnames(groupdf2mic)],groupdf2mic)
  gg <- ggplot(groupdf2comb, aes(fill=condition,y=Counts,x=condition)) +
    stat_summary(fun.y="mean", geom="bar") +
    facet_grid(Antibody~Gene, scales = "free_y") +
    geom_quasirandom(cex=0.5,varwidth = TRUE,dodge.width=0.3, groupOnX=TRUE) +
    labs(title = genelistname,x = "",y="normalized read counts") +
    #   scale_fill_manual(values=rev(brewer.pal(4,"RdYlBu"))) +
    scale_fill_manual(values=c("#0000FF","#008000")) +
    theme_bw() +
    theme(legend.direction = "vertical",legend.position="top",axis.text.x = element_text(angle = 45,vjust=0.5)) +
    scale_x_discrete(labels=gsub("\\+","\n\\+",levels(groupdf$condition))) 
  #   scale_y_continuous(breaks=seq(0, 4000, 2000), minor_breaks=seq(1000,3000,2000))
  ggsave(gg,file=paste0("read_counts_",genelistname,"_",label,".pdf"),width=1+0.5*length(genelist), height=5.3)
}
# ggsave(gg,file=paste0("Fig2_read_counts_orsgenes","_",label,".pdf"),width=3.5, height=5.3)
# ggsave(gg,file=paste0("Fig2_read_counts_orsgenes","_",label,".png"),width=3.5, height=5.3)

## finding DCB transcription factors
# stategy: look for pattern "transcription factor" in Notes field of gff file. 
# filter the ones that are DCB

save(stats_counts,mRNA,chip_comparisons,groupdfNoFlag,countsNoFlag,chr_sizes,genome,fls,gffgenes,febitA,febitF,chip_merged,chiplogFCsall,cnts,probenames,counts_edger,countsnall,groupdf,groupdf2,gwcov, cov.5e4,cov.5e5,cov.2e3,cov.5e4.norm,cov.2e3.norm,cov.5e5.norm,file="Robjects_nidulans.rdata")
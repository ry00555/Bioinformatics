#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
#options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

install.packages(optparse)
install.packages(data.table)
library(optparse)
library(data.table)

option_list = list(
    make_option(c("--sample_name"), dest="109_58", action="store"),
    make_option(c("--standardized_copy_ratios_file"), dest="/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.standardizedCR.tsv", action="store"),
    make_option(c("--denoised_copy_ratios_file), dest="/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.denoisedCR.tsv", action="store"),
    make_option(c("--contig_names"), dest="/home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict"[,2], action="store"),      #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths"), dest="/home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict"[,3]", action="store"),  #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--maximum_copy_ratio"), dest="maximum_copy_ratio", action="store", type="infinity"),
    make_option(c("--point_size_copy_ratio"), dest="point_size_copy_ratio", action="store", type="double"),
    make_option(c("--output_dir"), dest="/Users/rochelleyap/Desktop/LewisLab/Images", action="store"),
    make_option(c("--output_prefix"), dest="109_", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["109_58"]]
standardized_copy_ratios_file = opt[["/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.standardizedCR.tsv"]]
denoised_copy_ratios_file = opt[["/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.denoisedCR.tsv"]]
contig_names = opt[["/home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict"[,2]"]]
contig_lengths = opt[["/home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict"[,3]"]]
maximum_copy_ratio = opt[["maximum_copy_ratio"]]
point_size_copy_ratio = opt[["point_size_copy_ratio"]]
output_dir = opt[["/Users/rochelleyap/Desktop/LewisLab/Images"]]
output_prefix = opt[["109_"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.standardizedCR.tsv, /home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.denoisedCR.tsv, /home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict", /home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2.dict")))) {
    print(status=1))
}

GCF_000182925.2.dict"[,2] = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
GCF_000182925.2.dict"[,3] = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(GCF_000182925.2.dict"[,3])
contig_starts = c(0, head(contig_ends, -1))

CalculateMedianAbsoluteDeviation = function(dat) {
    return(median(abs(diff(dat))))
}

#plotting is extracted to a function for debugging purposes
WriteDenoisingPlots = function(109_58, /home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.standardizedCR.tsv, /home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.denoisedCR.tsv, GCF_000182925.2.dict"[,2], Run109CNV, 109_) {
    standardized_copy_ratios_df = ReadTSV(/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.standardizedCR.tsv)
    denoised_copy_ratios_df = ReadTSV(/home/ry00555/Bioinformatics/CrassaGenome/CopyRatios/109_58.denoisedCR.tsv)

    #transform to linear copy ratio
    standardized_copy_ratios_df[["COPY_RATIO"]] = 2^standardized_copy_ratios_df[["LOG2_COPY_RATIO"]]
    denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

    #determine copy-ratio midpoints
    standardized_copy_ratios_df[["MIDDLE"]] = round((standardized_copy_ratios_df[["START"]] + standardized_copy_ratios_df[["END"]]) / 2)
    denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

    #write the MAD files
    standardizedMAD = CalculateMedianAbsoluteDeviation(standardized_copy_ratios_df[["COPY_RATIO"]])
    denoisedMAD = CalculateMedianAbsoluteDeviation(denoised_copy_ratios_df[["COPY_RATIO"]])
    write.table(round(standardizedMAD, 3), file.path("/Users/rochelleyap/Desktop/LewisLab/Images", paste(Run109_, ".standardizedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(denoisedMAD, 3), file.path("/Users/rochelleyap/Desktop/LewisLab/Images", paste(Run109_, ".denoisedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(standardizedMAD - denoisedMAD, 3), file.path("/Users/rochelleyap/Desktop/LewisLab/Images", paste(Run109_, ".deltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((standardizedMAD - denoisedMAD) / standardizedMAD, 3), file.path("/Users/rochelleyap/Desktop/LewisLab/Images", paste(Run109_, ".scaledDeltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot standardized and denoised copy ratio on top of each other
    pre_color_blue = "#3B5DFF"
    post_color_green = "#4FC601"

    #plot up to maximum_copy_ratio (or full range, if maximum_copy_ratio = Infinity)
    denoising_plot_file = file.path("/Users/rochelleyap/Desktop/LewisLab/Images", paste(Run109_, ".denoised.png", sep=""))
    png(denoising_plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2, 1), cex=0.75, las=1)
    maximum_standardized_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(standardized_copy_ratios_df[["COPY_RATIO"]])
    SetUpPlot(109_58, "standardized copy ratio", 0, maximum_standardized_copy_ratio, paste("median absolute deviation = ", round(standardizedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatios(standardized_copy_ratios_df, pre_color_blue, contig_names, contig_starts, point_size_copy_ratio)
    maximum_denoised_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(denoised_copy_ratios_df[["COPY_RATIO"]])
    SetUpPlot(109_58, "denoised copy ratio", 0, maximum_denoised_copy_ratio, paste("median absolute deviation = ", round(denoisedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatios(denoised_copy_ratios_df, post_color_green, contig_names, contig_starts, point_size_copy_ratio)
    dev.off()

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(denoising_plot_file)))) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

#WriteDenoisingPlots(sample_name, standardized_copy_ratios_file, denoised_copy_ratios_file, contig_names, output_dir, output_prefix)

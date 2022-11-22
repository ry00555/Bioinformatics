#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=GATK.%j.out
#SBATCH --error=GATK.%j.err

cd $SLURM_SUBMIT_DIR
ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
ml picard/2.27.4-Java-13.0.2

#scp /Users/ry00555/Desktop/NeurosporaGenome/crassa.bed ry00555@@xfer.gacrc.uga.edu:/home/ry00555/Bioinformatics/CrassaGenome
#scp /Users/ry00555/Desktop/NeurosporaGenome/GCA_000182925.2_NC12_genomic.fasta ry00555@@xfer.gacrc.uga.edu:/home/ry00555/Bioinformatics/CrassaGenome
#curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > /home/ry00555/Bioinformatics/CrassaGenome/GCF_000182925.2_NC12_genomic.fna
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load BEDOPS/2.4.39-foss-2019b
ml Bowtie2/2.4.1-GCC-8.3.0
ml R/3.6.2-foss-2019b
#Converts the GFF file to BED format
samtools faidx GCF_000182925.2.fasta
bowtie2-build GCF_000182925.2.fasta reference
bowtie2 -x  reference -U GCF_000182925.2.fasta -S Crassa.sam

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  -R GCF_000182925.2.fasta \
  -O GCF_000182925.2.dict


java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
-I crassa.bed \
-R GCF_000182925.2.fasta \
-SD GCF_000182925.2.dict \
-O Crassa.interval_list

gatk PreprocessIntervals \
          -R GCF_000182925.2.fasta \
          -L Crassa.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          --bin-length 1000 \
          --padding 0 \
          -O Crassa.preprocessed_intervals.interval_list

java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
        I=/scratch/ry00555/OutputRun109/SortedBamFiles/109_58_Genomic.bam \
        MODE=SUMMARY

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      -I /scratch/ry00555/OutputRun109/SortedBamFiles/109_58_Genomic.bam \
      -O /home/ry00555/Bioinformatics/CrassaGenome/Run109ReadGroups/109_58_Genomic.bamoutput.bam \
      -RGID 1 \
      -RGLB lib1 \
      -RGPL illumina \
      -RGPU S34 \
      -RGSM 109_58
mkdir Run109ReadGroups
samtools index /home/ry00555/Bioinformatics/CrassaGenome/Run109ReadGroups/109_58_Genomic.bamoutput.bam

gatk CollectReadCounts \
     -I /home/ry00555/Bioinformatics/CrassaGenome/Run109ReadGroups/109_58_Genomic.bamoutput.bam \
     -L Crassa.interval_list \
     --interval-merging-rule OVERLAPPING_ONLY \
     -O /home/ry00555/Bioinformatics/CrassaGenome/109tsv/109_58.counts.tsv

gatk AnnotateIntervals \
          -R GCF_000182925.2.fasta \
          -L Crassa.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          -O Crassa_annotated_intervals.tsv

mkdir PanelofNormals
gatk CreateReadCountPanelOfNormals \
  -I 109tsv/109_58.counts.tsv \
  --annotated-intervals Crassa_annotated_intervals.tsv \
  -O PanelofNormals/109_58cnv.pon.hdf5

mkdir CopyRatios
  gatk DenoiseReadCounts \
          -I 109tsv/109_58.counts.tsv \
          --annotated-intervals Crassa_annotated_intervals.tsv \
          --count-panel-of-normals PanelofNormals/109_58cnv.pon.hdf5 \
          --standardized-copy-ratios CopyRatios/109_58.standardizedCR.tsv \
          --denoised-copy-ratios CopyRatios/109_58.denoisedCR.tsv

gatk PlotDenoisedCopyRatios \
                --standardized-copy-ratios CrassaGenome/CopyRatios/109_58.standardizedCR.tsv \
                --denoised-copy-ratios CrassaGenome/CopyRatios/109_58.denoisedCR.tsv \
                --sequence-dictionary CrassaGenome/GCF_000182925.2.dict \
                --output-prefix Run109CNV \
                --output CrassaGenome


scp -r ry00555@xfer.gacrc.uga.edu:/home/ry00555/Bioinformatics/CrassaGenome/Run109CNV.denoised.png ~/Downloads

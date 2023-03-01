#!/bin/bash
#SBATCH --job-name=j_CNVator
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=4
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=CNVnator.%j.out
#SBATCH --error=CNVnator.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )


##make output directory
OUTDIR= "/scratch/ry00555/Bioinformatics/CNVator"

module spider CNVnator/0.4.1-foss-2019b-ROOT-6.14.06

FILE=$1
SbPTH='/Users/ry00555/Desktop/RochelleLabDesktopData/IGV/mus30xmei3/Bamfiles'
REF='/Users/ry00555/Desktop/NeurosporaGenome/GCF_000182925.2.fasta'
CHROM='/Users/ry00555/Desktop/NeurosporaGenome/Chromsomes/*.fa'

#Stop at any error

cd /scratch/ry00555/Bioinformatics/CNVator
set -ue
#EXTRACTING READ MAPPING FROM BAM/SAM FILES
cnvnator -root ${FILE}.root -chrom $CHROM -tree ${SbPTH}/${FILE}.bam

  cnvnator -root ${FILE}.root -chrom $CHROM -his 1670 -fasta ${REF}
  cnvnator -root ${FILE}.root -chrom $CHROM -stat 1670
  cnvnator -root ${FILE}.root -chrom $CHROM -eval 1670
  cnvnator -root ${FILE}.root -chrom $CHROM -partition 1670
 #cnvnator -root ${FILE}.root -genotype 1670
 #cnvnator -root ${FILE}.root -view 1670
 cnvnator -root ${FILE}.root -chrom $CHROM -call 1670 > ${FILE}_cnv_longbin

 #First time I ran CNVnator I ran a smaller bin size (75). It will be run a second time but now with a bigger bin size 1670 (average gene lenght)

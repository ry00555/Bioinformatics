#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=4
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=CNVnator.%j.out
#SBATCH --error=CNVnator.%j.err

d $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )


##make output directory
OUTDIR= "/scratch/ry00555/Bioinformatics/CNVator"
if [ ! -d $OUTDIR ]
then
mkdir -p $OUTDIR
fi
module spider CNVnator/0.4.1-foss-2019b-ROOT-6.14.06

FILE=$1
SbPTH='/Users/ry00555/Desktop/RochelleLabDesktopData/IGV/mus30xmei3/mus30Samples'
REF='/Users/ry00555/Desktop/NeurosporaGenome/GCA_000182925.2_NC12_genomic.fasta'
#Stop at any error
set -ue
#move to CNVnator directory
cd /scratch/ry00555/Bioinformatics/CNVator
#EXTRACTING READ MAPPING FROM BAM/SAM FILES
nvnator -root ${FILE}.root -tree ${SbPTH}/${FILE}_Genomic.bam

 cnvnator -root ${FILE}.root -his 1670 -fasta ${REF}
 cnvnator -root ${FILE}.root -stat 1670
 cnvnator -root ${FILE}.root -eval 1670
 cnvnator -root ${FILE}.root -partition 1670
 #cnvnator -root ${FILE}.root -genotype 1670
 #cnvnator -root ${FILE}.root -view 1670
 cnvnator -root ${FILE}.root -call 1670 > ${FILE}_cnv_longbin

 #First time I ran CNVnator I ran a smaller bin size (75). It will be run a second time but now with a bigger bin size 1670 (average gene lenght)

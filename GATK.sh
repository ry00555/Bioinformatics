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
ml GATK/4.2.5.0-GCCcore-8.3.0-Java-1.8.lua



gatk CollectReadCounts \
          -I /scratch/ry00555/OutputRun109/SortedBamFile/109_58_Genomic.bam \
          -L intervals.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          -O 109_58.counts.tsv

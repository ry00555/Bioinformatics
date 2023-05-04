#!/bin/bash
#SBATCH --job-name=SRRtoGENE	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ry00555/Bioinformatics/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=ry00555@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, SUBMITTED, END, FAIL, ALL)
#SBATCH --output=../SRRtoGENE.%j.out
#SBATCH --error=../SRRtoGENE.%j.err

# create a directory for the sample names
mkdir /Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/SampleNames/

File = '/Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/SRRfromMetaBatch.txt'

# create a new file to save the sample names
touch /Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/SampleNames/sample_names.txt

# loop through each SRR ID in the list
while read srr; do
  # extract the sample name using vdb-dump and append it to the sample_names.txt file
  vdb-dump -S -f tab -C SAMPLE ~/ncbi/public/sra/SRR/$srr/$srr.sra | awk '{print $2}' >> /Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/SampleNames/sample_names.txt
done < /Users/ry00555/Desktop/RochelleLabDesktopData/RNAseq/OLDmappedJGIcounts/SRRfromMetaBatch.txt

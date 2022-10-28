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

#Set directory
cd /scratch/ry00555/Bioinformatics

#download metadata for an individual SRR sample
File='/scratch/ry00555/Bioinformatics/JGIAllCountsSRRONLY.txt'
for i in $(cat $File)
do
wget -O  $i.txt "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR$i&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,sample_alias,sample_title&format=tsv&download=true&limit=0"
done


#export only the gene name
grep -v study $i.txt | awk '{print $NF}' $i.txt >> TESTallSRRtoGENE.txt

#transfer files to RNAseq folder in RochelleLabDesktopData
#scp -r /scratch/ry00555/Bioinformatics/allSRRtoGENE.txt  $HOME/Desktop/RochelleLabDesktopData/RNAseq


#provided by Casey
#for i in SRR8269612

#do
#wget -O $i.txt "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$i&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,sample_alias,sample_title&format=tsv&download=true&limit=0"
#done

#grep -v study SRR*.txt
#awk '{print $NF}' SRR*.txt >> allSRRtoGENE.txt

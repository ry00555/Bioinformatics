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

cd /scratch/ry00555/Bioinformatics/

# Path to the file containing the list of SRR IDs
SRR_FILE="SRRfromMetaBatch.txt"

# Path to the file to write the sample names to
SAMPLE_NAMES_FILE="sample_names.txt"

# Loop through the SRR IDs in the file
while read SRR_ID; do
  # Download the HTML content of the SRA website for the SRR ID
  HTML=$(wget -qO- "https://www.ncbi.nlm.nih.gov/sra/?term=$SRR_ID")

  # Extract the sample name from the HTML content using grep
  SAMPLE_NAME=$(echo "$HTML" | grep -oP '(?<=<th scope="row">Sample:</th><td>)[^<]+' | head -1)

  # Print the SRR ID and sample name to the console
  echo "$SRR_ID, $SAMPLE_NAME"

  # Append the SRR ID and sample name to the sample names file
  echo "$SRR_ID, $SAMPLE_NAME" >> "sample_names.txt"
done < "$SRR_FILE"

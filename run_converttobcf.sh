#!/bin/bash

#SBATCH --account=def-lukens
#SBATCH --time=0-03:00:00
#SBATCH --ntasks-per-node=1

# Load required modules
module load bcftools

# Set paths to input and output files
fasta_file=../burbot_2021.fasta
bam_dir=$(pwd)
output_dir=$(pwd)

# Loop through all bam files in the input directory and generate a bcf file for each
for bam_file in ${bam_dir}/*.sorted.bam; do
  # Extract sample name from bam file name
  sample_name=$(basename ${bam_file} .bam)
  
  # Align the bam file to the reference fasta and generate a bcf file
  bcftools mpileup -f ${fasta_file} ${bam_file} -o ${output_dir}/${sample_name}.bcf
done


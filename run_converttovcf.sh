#!/bin/bash 

#SBATCH --account=def-lukens
#SBATCH --time=0-03:00:00
#SBATCH --ntasks-per-node=1

# Load required modules
module load bcftools

# Set paths
fasta_file=../burbot_2021.fasta 
bam_dir=$(pwd)
output_file=./filtered_bcftools.vcf

for bam_file in ${bam_dir}/*.sorted.bam
do
    sample_name=$(basename ${bam_file} .sorted.bam) 
    bcftools mpileup -d 50000 -f ${fasta_file} ${bam_file} | bcftools call -m --variants-only -o ${output_file}

done


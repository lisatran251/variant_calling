#!/bin/sh 

## This script convert .sam files to .bam files 
## Usage: sbatch run_bam.sh genome.sam 
## Note, will match .gz.sam 

#SBATCH --account=def-lukens 
#SBATCH --time=0-03:00:00 
#SBATCH --ntasks-per-node=1 

## Load modules 
module  load samtools 


## Converting .gz.sam files to .bam files 
for file in  *.gz.sam 
do 
	echo "Converting $file: " 
	## Remove '.gz.sam' extension and convert files to .bam format
	samtools view -S -b $file > ${file%.gz.sam}.bam 
done

echo "The process is now complete"

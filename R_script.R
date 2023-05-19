#install.packages("magrittr")
library(magrittr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

#####################################

## Allele Frequency 

# bcftools
var_freq_bcf <- read_delim("freq_bcftools.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# extract the first frequency column 
bcf_data <- var_freq_bcf$a1

# define a regular expression pattern to match numbers
pattern <- "\\d+\\.?\\d*"

bcf_freqs <- as.numeric(unlist(str_extract_all(bcf_data, pattern)))

# convert the extracted numbers to numeric format
bcf_freqs <- as.numeric(bcf_freqs)

# print the extracted numbers
print(bcf_freqs)

# plot the allele frequency 
hist(bcf_freqs, xlim=c(0,1), xlab= 'MAF', ylab='Number of SNPs', main='MAF Histogram of Bcftools')

#####################################

# freebayes
var_freq_freebayes <- read_delim("freq_freebayes.frq", delim = "\t",
                           col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# extract the first frequency column 
freebayes_data <- var_freq_freebayes$a1

# Define a regular expression pattern to match numbers
pattern <- "\\d+\\.?\\d*"

freebayes_freqs <- as.numeric(unlist(str_extract_all(freebayes_data, pattern)))

# Convert the extracted numbers to numeric format
freebayes_freqs <- as.numeric(freebayes_freqs)

# Print the extracted numbers
print(freebayes_freqs)

# Plot the allele frequency 
hist(freebayes_freqs, xlim=c(0,1), xlab= 'MAF', ylab='Number of SNPs', main='MAF Histogram of Freebayes')

#####################################

# bcftools 

# read the .txt file
depth_bcftools_file <- readLines("depth_perlocus_bcf.txt")

# extract the numbers after "DP="
depth_bcftools <- as.numeric(gsub("DP=", "", depth_bcftools_file))

# plot a histogram
hist(depth_bcftools,breaks = seq(0,3000, by=50), xlab='Mean depth per SNP', ylab='Number of SNPs', main='Depth Histogram of Bcftools')

#####################################

##freebayes

# read the .txt file
depth_freebayes_file <- readLines("depth_perlocus_freebayes.txt")

# extract the numbers after "DP="
depth_freebayes <- as.numeric(gsub("DP=", "", depth_freebayes_file))

# plot a histogram
hist(depth_freebayes,breaks = seq(0,6000, by=50), xlab='Mean depth per SNP', ylab='Number of SNPs', main='Depth Histogram of Freebayes')


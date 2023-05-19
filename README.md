# Variant Calling Genomic Project 
## Introduction
Variant calling is an essential step in next-generation sequencing data analysis. It helps researchers discover genomic variations, such as single nucleotide polymorphisms (SNPs) and insertions/deletions (Indels), between a sample and a reference genome. These variations provide valuable insights into disease genetics, population genetics, and evolutionary studies.

In this project, we compare two commonly used variant calling algorithms, freebayes and bcftools (or samtools). While both techniques employ probabilistic models to detect variations, they differ in their underlying algorithms and assumptions, resulting in tradeoffs. freebayes utilizes a Bayesian approach, considering the uncertainty of sequencing data and genotype calling, making it more accurate for rare variant detection in low-coverage regions. On the other hand, bcftools (or samtools) uses a maximum likelihood estimation approach, relying on a reference genome, which makes it more resource-efficient and suitable for handling large datasets with extensive coverage. However, it may be less accurate in low-coverage regions or miss rare mutations.

The objective of this study is to evaluate the performance of freebayes and bcftools in calling variants in ten individual burbot genomes with a burbot reference genome. By comparing the results, we aim to highlight the advantages and disadvantages of these algorithms.

## Method
For this study, we focus on single nucleotide polymorphisms (SNPs) obtained from ten individual burbot genomes and their corresponding reference genome. We use freebayes and bcftools as the variant calling methods.

Freebayes is chosen for its accuracy in detecting rare variants in low-coverage regions and its demonstrated high accuracy and sensitivity for SNP detection. Bcftools (or samtools) is selected for its efficiency in handling large datasets with extensive coverage and its previous success in SNP detection.

To ensure the accuracy of the SNP calling results, we apply strict filtering criteria. Variants with low-quality scores, low allele frequencies, and absence in at least 80% of the samples are excluded. This stringent filtering strategy aims to reduce false positive SNP calls and maintain the integrity of the SNP dataset. The combination of freebayes and bcftools, along with the filtering criteria, allows us to identify and evaluate SNPs in the burbot genomes with high accuracy.

The directory structure for this project is as follows:
Directory: /scratch/thuyduye/genomic_method_w23/genome_project3
Workflow details can be found at: /scratch/thuyduye/genomic_method_w23/genome_project3/workflow

## References
Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing.
Hwang, S., Kim, E., Lee, I., & Marcotte, E. M. (2015). Systematic comparison of variant calling pipelines using gold standard personal exome variants. Scientific Reports, 5(1), 17875. https://doi.org/10.1038/srep17875
Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

## Author
This project was developed by Lisa Tran.

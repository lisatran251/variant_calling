# run .sh script to convert sam to bam files 
sbatch convertsamtobam.sh

# sort .bam files 
find . -type f -name '*.bam' -exec sh -c 'samtools sort -o "${1%.bam}.sorted.bam" "$1"' _ {} \;

# convert sorted.bam to bcf file
bcftools mpileup -a DP,AD --skip-indels -P ILLUMINA -f ../burbot_2021.fasta *.sorted.bam -o burbot_bcftools.bcf 
##wating to load, taking forever 

#convert bcf file to vcf file 
bcftools call -m --variants-only --format-fields GQ --skip-variants indels burbot_bcftools.bcf | bcftools filter --set-GTs . --include 'QUAL > 20' |  bcftools view --min-alleles 2 --max-alleles 2 --types snps --apply-filter "PASS" --output-type v --output-file burbot_bcftools_filtered.vcf
## at this point only one vcf file contains all samples 

#add read group 
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_70_TRL_363q.sorted.bam O=871_70_TRL_363q_rg.sorted.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small1
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_71_TRL_421q.sorted.bam  O=871_71_TRL_421q_rg.sorted.bam  RGID=2 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small2
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_72_TRL_571q.sorted.bam   O=871_72_TRL_571q_rg.sorted.bam   RGID=3 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small3
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_73_TRL_571q.sorted.bam   O=871_73_TRL_571q_rg.sorted.bam   RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small4
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_74_TRL_599q.sorted.bam   O=871_74_TRL_599q_rg.sorted.bam   RGID=5 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small5
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_75_TRL_449q.sorted.bam   O=871_75_TRL_449q_rg.sorted.bam   RGID=6 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small6
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_76_TOL_358q.sorted.bam   O=871_76_TOL_358q_rg.sorted.bam   RGID=7 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small7
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_77_TOL_360q.sorted.bam   O=871_77_TOL_360q_rg.sorted.bam   RGID=8 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small8
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_78_TOL_325q.sorted.bam   O=871_78_TOL_325q_rg.sorted.bam   RGID=9 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small9
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=871_79_TOL_355q.sorted.bam   O=871_79_TOL_355q_rg.sorted.bam   RGID=10 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=small10


#indexing files 
samtools index 871_70_TRL_363q_rg.sorted.bam
samtools index 871_71_TRL_421q_rg.sorted.bam
samtools index 871_72_TRL_571q_rg.sorted.bam
samtools index 871_73_TRL_571q_rg.sorted.bam
samtools index 871_74_TRL_599q_rg.sorted.bam
samtools index 871_75_TRL_449q_rg.sorted.bam
samtools index 871_76_TOL_358q_rg.sorted.bam
samtools index 871_77_TOL_360q_rg.sorted.bam
samtools index 871_78_TOL_325q_rg.sorted.bam
samtools index 871_79_TOL_355q_rg.sorted.bam

#freebayes 
freebayes -f ../burbot_2021.fasta -b 871_70_TRL_363q_rg.sorted.bam -b 871_71_TRL_421q_rg.sorted.bam -b 871_72_TRL_571q_rg.sorted.bam -b 871_73_TRL_571q_rg.sorted.bam -b 871_74_TRL_599q_rg.sorted.bam -b 871_75_TRL_449q_rg.sorted.bam -b 871_76_TOL_358q_rg.sorted.bam -b 871_77_TOL_360q_rg.sorted.bam -b 871_78_TOL_325q_rg.sorted.bam -b 871_79_TOL_355q_rg.sorted.bam | vcffilter -f "QUAL > 20"  > burbot_bcftools_simple_2.vcf

#simple_bcftools 
bcftools mpileup -f ../burbot_2021.fasta *sorted.bam | bcftools call -m --variants-only > burbot_bcftools_simple.vcf

#add filter for freebayes (additional) 
freebayes -f ../burbot_2021.fasta -b 871_70_TRL_363q.bam -b 871_71_TRL_421q.bam -b 871_72_TRL_571q.bam -b 871_73_TRL_571q.bam -b 871_74_TRL_599q.bam -b 871_75_TRL_449q.bam -b 871_76_TOL_358q.bam -b 871_77_TOL_360q.bam -b 871_78_TOL_325q.bam -b 871_79_TOL_355q.bam | vcffilter -f "QUAL>20" > ex.vcf


#Run vcftools to compare 2 vcfs #This one is with QUAL > 20 only 
vcftools --vcf burbot_freebayes.vcf --diff burbot_bcftools_filtered.vcf --out vcf_compare_vcftools


## apply --mac filter to remove monomorphic SNPs. Since monomorphic SNPs have allele frequencies of 0 and 1, —mac 2 filter will exclude variants with a minor allele count less than 2
vcftools --vcf burbot_freebayes.vcf --mac 2  --recode  --recode-INFO-all --out burbot_freebayes_mac_2

vcftools --vcf burbot_bcftools_filtered.vcf --mac 2  --recode  --recode-INFO-all --out burbot_bcftools_mac_2

## Use the new vcf files (burbot_freebayes_mac_2.recode.vcf and burbot_bcftools_mac_2.recode.vcf) for other comparisions 
bcftools stats burbot_bcftools_mac_2.recode.vcf  | grep "number of SNPs"
bcftools stats burbot_freebayes_mac_2.recode.vcf  | grep "number of SNPs"

## Comparing overlap between VCFs
vcftools --vcf burbot_freebayes_mac_2.recode.vcf --diff burbot_bcftools_mac_2.recode.vcf --diff-site --not-chr GWHBAUG00000009 --not-chr GWHBAUG00000008 --out vcf_compare_vcftoolsvcftools --vcf burbot_freebayes_mac_2.recode.vcf --diff burbot_bcftools_mac_2.recode.vcf --diff-site --not-chr GWHBAUG00000009 --not-chr GWHBAUG00000008 --out vcf_compare_vcftools

## Extract the per-locus "DP" field 
grep -oP 'DP=\d+' burbot_bcftools_mac_2.recode.vcf > depth_perlocus_bcf.txt
grep -oP 'DP=\d+' burbot_freebayes_mac_2.recode.vcf > depth_perlocus_freebayes.txt

## Compare overlap between vcf 
bgzip burbot_freebayes_mac_2.recode.vcf
bgzip burbot_bcftools_mac_2.recode.vcf 

bcftools index burbot_freebayes_mac_2.recode.vcf.gz
bcftools index burbot_bcftools_mac_2.recode.vcf.gz

bcftools isec -p isec burbot_bcftools_mac_2.recode.vcf.gz burbot_freebayes_mac_2.recode.vcf.gz

#Calculate allele depth 
vcftools --vcf burbot_bcftools_mac_2.recode.vcf --site-depth > alleledept_locus_bcftools 
vcftools --vcf burbot_freebayes_mac_2.recode.vcf --site-depth > alleledept_locus_freebayes

#Calculate individual depth 
vcftools --vcf burbot_bcftools_mac_2.recode.vcf --depth -c > ind_depth_bcftools.txt
vcftools --vcf burbot_freebayes_mac_2.recode.vcf --depth -c > ind_depth_freebayes.txt

#Calculate allele frequency for ea locus
vcftools --vcf burbot_bcftools_mac_2.recode.vcf --freq --out freq_bcftools 
vcftools --vcf burbot_freebayes_mac_2.recode.vcf --freq --out freq_freebayes


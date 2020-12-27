#!/bin/bash
#
# Downloading reference files
wget -nc  http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit
wget -nc  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
./twoBitToFa hg19.2bit hg19.fa
mv hg19.fa genome.fa
samtools faidx genome.fa
bwa index genome.fa
gatk CreateSequenceDictionary -R genome.fa

# Downloading dbSNP files
wget -nc ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
gunzip -c All_20180423.vcf.gz > All_20180423.vcf
cat All_20180423.vcf | parallel --pipe sed '/^#/!s/^/chr/' > dbSNP_hg19.vcf
rm All_20180423.vcf
rm All_20180423.vcf.gz
gatk IndexFeatureFile -F dbSNP_hg19.vcf

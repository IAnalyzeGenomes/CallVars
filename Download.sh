#!/bin/bash
#
set -e
# Downloading human reference genome files 
wget -nc  http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit
wget -nc  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
./twoBitToFa hg19.2bit hg19.fa
mv hg19.fa genome.fa

# Getting human reference genome files ready for analysis
samtools faidx genome.fa
bwa index -b 50000000 genome.fa
gatk --java-options -Xmx4g CreateSequenceDictionary -R genome.fa

# Downloading and getting dbSNP files ready 
wget -nc ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
zcat All_20180423.vcf.gz | parallel --pipe sed '/^#/!s/^/chr/' > dbSNP_hg19.vcf
rm All_20180423.vcf.gz
gatk --java-options -Xmx4g SortVcf -I dbSNP_hg19.vcf -O dbSNP_hg19_sorted.vcf  -SD genome.dict
rm dbSNP_hg19.vcf

#!/bin/bash
#
# Downloading human reference genome files 
wget -nc  http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit
wget -nc  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
./twoBitToFa hg19.2bit hg19.fa
mv hg19.fa genome.fa

## Getting human reference genome files ready for analysis
samtools faidx genome.fa
RESULT1=$?
if [ $RESULT1 -eq 0 ]; then
  echo samtools_faidx_success
else
  echo samtools_faidx_failed
fi
bwa index genome.fa
RESULT2=$?
if [ $RESULT2 -eq 0 ]; then
  echo bwa_index_success
else
  echo bwa_index_failed
fi
gatk CreateSequenceDictionary -R genome.fa
RESULT3=$?
if [ $RESULT3 -eq 0 ]; then
  echo gatk_CreateSequenceDictionary_success
else
  echo gatk_CreateSequenceDictionary_failed
fi

# Downloading and getting dbSNP files ready 
wget -nc ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
zcat All_20180423.vcf.gz | parallel --pipe sed '/^#/!s/^/chr/' > dbSNP_hg19.vcf
RESULT4=$?
if [ $RESULT4 -eq 0 ]; then
  echo dbSNP_Download_success
else
  echo dbSNP_Download_failed
fi
rm All_20180423.vcf.gz
gatk IndexFeatureFile -F dbSNP_hg19.vcf
RESULT5=$?
if [ $RESULT5 -eq 0 ]; then
  echo gatk_IndexFeatureFile_success
else
  echo gatk_IndexFeatureFile_failed
fi

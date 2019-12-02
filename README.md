# CallVars: 

CallVars is an automated, reproducible Snakemake workflow taking Illumina paired-end FastQ files directly to a filtered list of high confidence variants for clinical review. This workflow largely follows Broad Institute's "Best Practices" guidlines for germline short variant discovery (SNPs + Indels) for single sample and also reports a filtered list of somatic variants. 

CallVars can be helpful to anyone working with targeted gene panels or whole exomes for rare disease or cancer diagnosis/treatment. If you think CallVars can help with your study feel free to DM me on twitter (@IAnalyzeGenomes) with questions. Also, feedback, comments or bug reports are welcome for me to integrate improvements.

CallVars currently uses GnomAD allele frequency as a key filter to report variants, having either genomes or exomes allele frequency less than 0.5%, for clinical review. Also, below GATK guidelines were used to apply generic hard-filtering to PASS/FAIL a variant.  
https://software.broadinstitute.org/gatk/documentation/article.php?id=6925

# CallVars workflow overview:
CallVars sequentially performs below steps of Next-Gen Sequencing (NGS) analysis.
1) Pre-processing using Cutadapt
2) Mapping using BWA
3) Sorting using samtools
4) Removing duplicates using GATK MarkDuplicates
5) Indexing using samtools
6) Base quality score recalibration using GATK BaseRecalibrator and ApplyBQSR
7) Germline variant detection using GATK HaplotypeCaller
8) Somatic variant detection using GATK Mutect2
9) Functional annotation for germline variants using GATK Funcotator
10) Functional annotation for somatic variants using GATK Funcotator
11) Variant filtration for germline variants using GATK VariantFiltration
12) Variant filtration for somatic variants using GATK VariantFiltration

# Setting up a working directory to run CallVars:
All the below listed files/folders must be present in the working directory before you run CallVars. Make sure the names of files/folders match exctly as listed.  

	- "FastQ" folder containing paired-end reads ending in _R1.fastq and _R2.fastq (see attached for test files A_R1.fastq and A_R2.fastq)
	- "CallVars.yml" (attached in repo)
	- "Snakefile" (attached in repo)
	- "Target.bed" (Your target file in BED format)

You will need to download each of the below listed files from their respective public repositories. I have provided links to these resources.

	- "1000G_phase1.indels.hg19.sites.vcf"
	- "1000G_phase1.indels.hg19.sites.vcf.idx"
	- "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
	- "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
	These 1000 genome indel files can be downloaded using below link.
	ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
	
	- "dbSNP_20180423.vcf"
	- "dbSNP_20180423.vcf.idx"
	These dbSNP files can be downloaded usibg below link.
	ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
	
	- "GNOMAD_hg19.vcf"
	- "GNOMAD_hg19.vcf.idx"
	These gnomAD files can be downloaded using below link.
	http://hgdownload.cse.ucsc.edu/gbdb/hg19/gnomAD/vcf/
	
	- HG19 Reference Genome Folder ‘UCSCWholeGenomeFasta’ containing files
  	"genome.dict"
	"genome.fa"
	"genome.fa.amb" 
	"genome.fa.ann"
	"genome.fa.bwt"
	"genome.fa.fai"
	"genome.fa.pac"
	"genome.fa.sa"
	"GenomeSize.xml"
	The reference genome files can be downloaded in 2bit format using below link.
	http://hgdownload.cse.ucsc.edu/gbdb/hg19/
	The utility program, twoBitToFa (available from the kent src tree), can be used to extract .fa file(s) from this file.  A pre-compiled version of the command line tool can be
    found at:
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
	Converting 2bit to fa --> ./twoBitToFa hg19.2bit hg19.fa
	renaming --> mv hg19.fa genome.fa
	indexing .fa file --> bwa index genome.fa
	Creating .dict file --> gatk CreateSequenceDictionary –R genome.fa

	
	- dataSourcesFolder containing below data source. 
	Genecode, Clinvar, Gnomad
	Data source was downloaded using below link.
	https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124s.tar.gz

# Running CallVars on Linux terminal: 
1)	Check the working directory and FastQ files: 

In the linux terminal, change the directory to the working directory that contains all the needed files and folders for running snakemake.
Make sure your fastq files are in a ‘FastQ’ directory and they end in ‘_R1.fastq’ and ‘_R2.fastq’, say A_R1.fastq and A_R2.fastq. The pipeline also works with gzipped fastq files, say A_R1.fastq.gz and A_R2.fastq.gz.

2)	Install miniconda: 

Use below link to install miniconda.

https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

3)	Install snakemake:
	
		conda install -c bioconda -c conda-forge snakemake

4)	Create environment CallVars:
	
		conda env create –n CallVars –f CallVars.yml

5)	Activate CallVars environment:

		conda activate CallVars

6)	Running snakemake: 
		Ensure you run the below command in working directory.

		snakemake CallVars_Output/Results/{your_sample_name}_CallVars_Germline.txt CallVars_Output/VCF/{your_sample_name}_CallVars_Somatic.vcf --cores N

For instance, if the names if the FastQ files are A_R1.fastq and A_R2.fastq and number of cores available are 6 then run below command.

		snakemake CallVars_Output/Results/A_CallVars_Germline.txt CallVars_Output/VCF/A_CallVars_Somatic.txt --cores 6

After the worklow has run successfully, below listed files will be available for clinician's review.

	1] CallVars_Output/Results/A_CallVars_Germline.txt containing a filtered list of germline variants.
	2] CallVars_Output/Results/A_CallVars_Somatic.txt containing a filtered list of somatic variants.
	3] CallVars_Output/VCF/A_germline_func_filter.vcf containing a complete list of germline variants.
	4] CallVars_Output/VCF/A_somatic_func_filter.vcf containing a complete list of somatic variants.


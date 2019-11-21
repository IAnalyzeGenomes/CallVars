# CallVars: 

CallVars is an automated reproducible Snakemake workflow that takes you from Illumina paired-end reads (FastQ files) to functionally annotated variants (VCF files) for human whole exome sequencing or targeted sequencing data. This workflow largely follows Broad Institute's "Best Practices" guidlines for germline short variant discovery (SNPs + Indels) in single sample and can also detect somatic variants. 

This workflow has been tested using 64-bit linux OS. Note that it uses hg19 version of human reference genome for analysis of Illumina paired-end reads. Below listed analysis steps are performed by the workflow sequentially. 

# Steps in Snakemake workflow:
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

# Setting up working directory to run Snakemake:
All the below listed files need to be present in the working directory before you can run the snakemake workflow.

	- ‘FastQ’ directory containing paired-end reads (ending in _R1.fastq and _R2.fastq, say A_R1.fastq and A_R2.fastq) 
	- CallVars.yml
	- Snakefile 

Note that the workflow currently uses hg19 version of human reference genome. You will need to download each of the below listed files from their respective online public repositories. However, if you are not able to find any of the below listed resources, feel free to reach me at amitbinf[at]med.umich.edu or amit4biotek[at]gmail.com for assistance.

	- 1000G_phase1.indels.hg19.sites.vcf
	- 1000G_phase1.indels.hg19.sites.vcf.idx
	- dbSNP_20180423.vcf
	- dbSNP_20180423.vcf.idx
	- GNOMAD_hg19.vcf
	- GNOMAD_hg19.vcf.idx
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
	- HG19 Reference Genome Directory ‘UCSCWholeGenomeFasta’ containing files
  	genome.dict, genome.fa, genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.fai, genome.fa.pac, genome.fa.sa, GenomeSize.xml
	- dataSourcesFolder containing below data source. 
	Genecode, Clinvar, Gnomad
	Data source was downloaded using below link.
	https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124s.tar.gz
# Running Snakemake workflow on Linux terminal: 
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

		snakemake CallVars_Output/VCF/{your_sample_name}_germline.vcf CallVars_Output/VCF/{your_sample_name}_somatic_func.vcf --cores N

For instance, if the names if the FastQ files are A_R1.fastq and A_R2.fastq and number of cores available are 6 then run below command.

		snakemake CallVars_Output/VCF/A_germline.vcf CallVars_Output/VCF/A_somatic_func.vcf --cores 6

Where CallVars_Output/VCF/A_germline.vcf is the germline variant file, CallVars_Output/VCF/A_somatic.vcf is the somatic variant file and N is number of cores.

 I plan to continue refining this workflow to detect manageable list of clinically relevant variants for whole exome sequencing or targeted sequencing data. Feel free to reach me at amitbinf[at]med.umich.edu or amit4biotek[at]gmail.com with any questions.

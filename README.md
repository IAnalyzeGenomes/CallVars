# CallVars: 

CallVars is an automated, reproducible Snakemake workflow that takes Illumina paired-end FastQ files directly to a handful of high confidence variants for clinical review. This workflow largely follows Broad Institute's "Best Practices" guidlines for germline short variant discovery (SNPs + Indels) for single sample and also reports a filtered list of somatic variants. CallVars can be helpful to anyone using targeted gene panels or whole exomes for rare disease or cancer diagnosis/treatment.

GnomAD allele frequency is a key filter used to report variants, having either genomes or exomes allele frequency less than 0.5%, for clinical review. Below GATK guidelines were used to apply generic hard-filtering to PASS/FAIL a variant.  
https://software.broadinstitute.org/gatk/documentation/article.php?id=6925

Note that CallVars uses hg19 version of human reference genome for Next-Gen Sequencing (NGS) data analysis. Below listed NGS analysis steps are performed by CallVars sequentially. 

# Steps/Rules in CallVars:
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


# Setting up working directory to run Snakemake:
Below I have provided step-by-step instructions on successfully running CallVars, however feel free to DM me on twitter (@IAnalyzeGenomes) if you need assistance in setting up and running this workflow.

All the below listed files/folders need must be present in the working directory before you can run CallVars. 

	- ‘FastQ’ folder containing paired-end reads ending in _R1.fastq and _R2.fastq (see attached for test files A_R1.fastq and A_R2.fastq)
	- CallVars.yml (see attached files)
	- Snakefile (see attached files)

You will need to download each of the below listed files from their respective online public repositories. I have provided links to these resources.

	- 1000G_phase1.indels.hg19.sites.vcf
	- 1000G_phase1.indels.hg19.sites.vcf.idx
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
	These files can be downloaded using below link.
	ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
	
	- dbSNP_20180423.vcf
	- dbSNP_20180423.vcf.idx
	These files can be downloaded usibg below link.
	ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
	
	- GNOMAD_hg19.vcf
	- GNOMAD_hg19.vcf.idx
	These files can be downloaded usibg below link.
	http://hgdownload.cse.ucsc.edu/gbdb/hg19/gnomAD/vcf/
	
	- HG19 Reference Genome Folder ‘UCSCWholeGenomeFasta’ containing files
  	genome.dict, genome.fa, genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.fai, genome.fa.pac, genome.fa.sa, GenomeSize.xml
	These files can be downloaded usibg below link.
	http://hgdownload.cse.ucsc.edu/gbdb/hg19/
	
	- dataSourcesFolder containing below data source. 
	Genecode, Clinvar, Gnomad
	Data source was downloaded using below link.
	https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124s.tar.gz

# Running Snakemake workflow using linux terminal: 
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

Where CallVars_Output/Results/A_CallVars_Germline.txt is the filtered germline variant file and CallVars_Output/VCF/A_CallVars_Somatic.txt is the filtered somatic variant file.

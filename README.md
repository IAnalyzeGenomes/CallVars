# CallVars:

**CallVars is an automated, reproducible and scalable Snakemake workflow that takes paired-end FastQ files directly to a filtered list of high confidence variants for clinical review. The workflow is based on [Broad Institute's Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) guidelines for germline short variant discovery (SNPs + Indels) for single sample.** 

CallVars can be helpful to anyone working with targeted gene panels or even whole exomes to detect disease causing/associated variants that can potentially help clinicians with a diagnosis/treatment. If CallVars can be of your help or if you run into issues while using it, feel free to DM me on twitter [(@IAnalyzeGenomes)](https://twitter.com/IAnalyzeGenomes). CallVars is under development so any feedback for its improvement is much appreciated.

CallVars is configured to run with parameters listed in "config.yaml" file, which is attached in this repository. You may change parameter values in config file to customize the workflow to your needs. You can list your samples (within "SAMPLE" section) in config file to scale the workflow. You can also choose to run the workflow either with hg19 or hg38 version of human reference genome. The description below is pertaining to hg19.

## Installing and running CallVars on linux Command Line Interface (CLI):

#### 1) Setting up a CallVars working directory:
		
		git pull https://github.com/IAnalyzeGenomes/CallVars

Pull this repo using above command and make sure all the below listed files/folders must be present in the working directory before you run CallVars.

	- "FastQ" folder containing paired-end reads ending in _R1.fastq and _R2.fastq (test files A_R1.fastq and A_R2.fastq)
	- "CallVars.yml" (This file is needed while creating a conda environment for running CallVars.)
	- "Snakefile" (This file is needed for running CallVars.)
	- "config.yaml" (This file lists samples and workflow parameter values. You may adjust these per your need.)
	- "Target.bed" (This file can be replaced by your target file of interest in BED format.)
	- "gatkPythonPackageArchive.zip" (This file is needed while creating a conda environment for running CallVars.)

**Please get the "dataSourcesFolder" directory ready. A google account will be required for download.**
	
	- "dataSourcesFolder" containing below data sources. 
			Genecode, Clinvar, Gnomad
	  These data sources can be downloaded using below link.
	  https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124g.tar.gz

#### 2) Install miniconda:

Use below link to install miniconda.

https://conda.io/projects/conda/en/latest/user-guide/install/linux.html


#### 3)	Setting up CallVars environment using conda:
	
		conda env create –n CallVars –f CallVars.yml

#### 4)	Activate CallVars environment using conda:

		conda activate CallVars

#### 5) Running the download script, Download.sh (attached in this repo)
		sh Download.sh

The above bash script can take 1.5-2 hours to run as it will download below files needed for CallVars workflow to run successfully.
		
		- dbSNP_hg19.vcf
		- dbSNP_hg19.vcf.idx
		- genome.dict
		- genome.fa
		- genome.fa.amb
		- genome.fa.ann
		- genome.fa.bwt
		- genome.fa.fai
		- genome.fa.pac
		- genome.fa.sa

#### 6)	Running CallVars: 
Ensure you run the below command's in the working directory.

Use below command on CLI for a dry run:
		
		snakemake -np

		
Use below command on CLI to execute the workflow:
		
		snakemake
		
Use below command on CLI if your machine supports multiple CPU cores:
		
		snakemake --cores N

**After the worklow has run successfully, below listed files will be available for clinical review. If you are running more samples, then you will see these files for all your samples.**

**1] CallVars/Reports/A_Germline.txt** contains a filtered list of germline variants from GATK HaplotypeCaller.

**2] CallVars/Reports/A_Germline_All.vcf** contains a complete list of germline variants from GATK HaplotypeCaller.

**3] CallVars/Reports/A_samtools.vcf** contains a full list variants from samtools mpileup.

**4] CallVars/NoDupReads/A_PerBaseCov.txt** contains coverage for each base.

**5] CallVars/NoDupReads/A_PerBaseCov_LessThan20.txt** contains bases for which coverage is less than 20.


# CallVars Workflow:
CallVars sequentially performs below steps of Next-Gen Sequencing (NGS) analysis.
![](CallVarsWorkflowEngine.PNG)


### 1) Pre-processing using Cutadapt

	rules/AdapterTrim.py

Pre-processing prepares the data for NGS analysis. When DNA or RNA molecules are sequenced using Illumina short reads technology, the machine may sequence into the adapter ligated to the 3’ end of each molecule during library preparation. Consequently, the reads that are output contain the sequence of the molecule of interest and also the adapter sequence. Also, with Illumina sequencing machines, the quality of reads is high at the beginning but degrades towards the 3’ end of the reads. 
	
CallVars uses Cutadapt to remove adapters from sequencing reads. Cutadapt also trims the read ends with quality below 20 and removes the ambiguous bases (N’s) from the reads ends. 

The adapters to be trimmed and threshold for quality values can be customized using the config.yaml file attached in this repository. 

### 2) Mapping using BWA

	rules/Mapping.py 

Once the high quality reads are obtained from pre-processing, the next step is mapping them to human reference genome. CallVars uses BWA-mem to map short Illumina paired-end reads to hg19/hg38 version of human reference genome. The reference genome to be used can be configured using the config.yaml file attached in the repo. 

This step generate a Binary Alignment Map also called a BAM file. The reference genome files needed for the analysis were downloaded in 2bit format using below link.
	
http://hgdownload.cse.ucsc.edu/gbdb/hg19/

### 3) BAM Preparation
	
	rules/BamPrep.py 
	
#### 	3A - Sorting using samtools
Now that we have a BAM file, we need to index it. All BAM files need an index, as they tend to be large and the index allows us to perform computationally complex operations on these files without it taking days to complete. Before we index the BAM file we need to sort them by position and remove duplicates. This step performs sorting the BAM file by position.
	
#### 	3B - Removing duplicates using GATK MarkDuplicates
CallVars uses this tool to locate and remove duplicate reads in a BAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. The MarkDuplicates tool works by comparing sequences in the 5 prime positions of read-pairs in a BAM file.
	
#### 	3C - Indexing using samtools
CallVars now performs indexing. 

### 4) 	Variant calling using samtools
	
	rules/SamtoolsVC.py

CallVars reports a list of variants from samtools mpileup and bcftools. Samtools collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format (see below). It does not call variants.Bcftools applies the prior and does the actual calling. It can also concatenate BCF files, index BCFs for fast random access and convert BCF to VCF. 

### 5) 	Coverage analysis using bedtools

	rules/PerBaseCov.py 

CallVars uses bedtools coverage utility to report per base coverage for the targeted file you provided in BED format. 
It also reports bases with coverage less than 20 to review regions of low coverage.

### 6) Germline variant detection using GATK HaplotypeCaller

	rules/HaplotypeCaller.py

CallVars uses GATK HaplotypeCaller to call germline SNPs and indels via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper. 

### 7) Functional annotation using GATK Funcotator

	rules/Funcotator.py

CallVars uses GATK Funcotator (FUNCtional annOTATOR) to analyze given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file. This tool is a functional annotation tool that allows a user to add annotations to called variants based on a set of data sources, each with its own matching criteria.
	
Data from Genecode, Clinvar and Gnomad were used for annotation of variants. Data source was downloaded using below link.
	
https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124g.tar.gz
	

### 8) Variant filtration using GATK VariantFiltration

	rules/VariantFiltration.py

[GATK guidelines](https://software.broadinstitute.org/gatk/documentation/article.php?id=6925) were used to apply generic hard-filtering to add PASS/FAIL tags to variants. Note that CallVars doesn't filter the variants based on PASS/FAIL tags. 
CallVars currently uses gnomAD allele frequency as a key filter to report variants having either genomes or exomes allele frequency less than 1% for clinical review. This value can be customized using the config.yaml file attached in this repository. 

## Benchmarking:
CallVars reported a 100% sensitivity for SNPs and 98.5% sensitivity for SNPs+Indels in four Genome In A Bottle (GIAB) samples [NA12878, NA24385, NA24143 and NA24149] combined, for a targeted panel of 64 cancer specific genes as listed below. VCFEVAL utility from Real Time Genomics was used to evaluate the sensitivity.
		
	ALK, APC, ATM, AXIN2, BAP1, BARD1, BMPR1A, BRCA1, BRCA2, BRIP1, CDC73, CDH1, CDK4, CDKN1C, CDKN2A, CHEK2,
	DICER1, EPCAM, FANCC, FH, FLCN, GPC3, GREM1, HOXB13, MAX, MEN1, MET, MITF, MLH1, MSH2, MSH6, MUTYH, NBN, 
	NF1, NF2, PALB2, PHOX2B, PMS1, PMS2, POLD1, POLE, PRKAR1A, PTCH1, PTEN, RAD51C, RAD51D, RB1, RET, SDHA, 
	SDHAF2, SDHB, SDHC, SDHD, SMAD4, SPRED1, STK11, SUFU, TMEM127, TP53, TSC1, TSC2, VHL, WT1 and XRCC2



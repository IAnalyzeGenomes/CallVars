# CallVars: 

CallVars is a Snakemake workflow that currently performs preliminary Next-Gen Sequencing (NGS) analysis i.e pre-processing, mapping and variant calling using Illumina paired-end reads (fastq files). The end goal of the workflow is to report a list of clinically relevant variants in targeted gene panels for human data which can be reveiewed by clinicians for disease diagnosis and/or selecting possible treatment options for the patient. 

The workflow uses hg19 version of human reference genome for analysis and it has been tested using 64-bit linux OS.

# Steps in Snakemake workflow:
1) Pre-processing using Cutadapt
2) Mapping using BWA
3) Sorting using samtools
4) Removing duplicates using GATK MarkDuplicates
5) Indexing using samtools
6) Base quality score recalibration using GATK BaseRecalibrator and ApplyBQSR
7) Germline variant detection using GATK HaplotypeCaller
8) Somatic variant detection using GATK Mutect2

# Setting up working directory to run Snakemake:
All the below listed files need to be present in the working directory before you can run the NGS analysis.

	- Paired-end reads (ending in _R1.fastq and _R2.fastq, say A_R1.fastq and A_R2.fastq) in ‘FastQ’ directory
	- CallVars.yml
	- Snakefile

Note that the workflow currently uses hg19 version of human reference genome. You will need to download each of the below listed files from their respective online public repositories. However, if you are not able to find any of the below listed files, feel free to reach me at amitbinf[at]med.umich.edu or amit4biotek[at]gmail.com.

	- 1000G_phase1.indels.hg19.sites.vcf
	- 1000G_phase1.indels.hg19.sites.vcf.idx
	- dbSNP_20180423.vcf
	- dbSNP_20180423.vcf.idx
	- GNOMAD_hg19.vcf
	- GNOMAD_hg19.vcf.idx
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	- Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
	- Reference Genome Directory ‘UCSCWholeGenomeFasta’ containing files
  	genome.dict, genome.fa, genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.fai, genome.fa.pac, genome.fa.sa, GenomeSize.xml

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

		snakemake CallVars_Output/VCF/{your_sample_name}_germline.vcf CallVars_Output/VCF/{your_sample_name}_somatic.vcf --cores N

For instance, if the names if the FastQ files are A_R1.fastq and A_R2.fastq and number of cores available are 6 then run below command.

		snakemake CallVars_Output/VCF/A_germline.vcf CallVars_Output/VCF/A_somatic.vcf --cores 6

Where CallVars_Output/VCF/A_germline.vcf is the germline variant file, CallVars_Output/VCF/A_somatic.vcf is the somatic variant file and N is number of cores.

 I plan to continue refining this workflow for detecting clinically relevant variants in targeted gene panels for human data. Feel free to reach me at amitbinf[at]med.umich.edu or amit4biotek[at]gmail.com with any questions.

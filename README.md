# CallVars
A snakemake variant calling pipeline for a single sample Illumina short paired-end reads


This is a snakemake pipeline for doing preliminary NGS analysis of converting raw Illumina short paired-end reads to a list of somatic and/or germline variants. 

Snakemake NGS analysis workflow:
1) Pre-processing using Cutadapt
2) Mapping using BWA
3) Sorting using samtools
4) Removing duplicates using GATK MarkDuplicates
5) Indexing using samtools
6) Base quality score recalibration using GATK BaseRecalibrator and ApplyBQSR
7) Germline variant detection using GATK HaplotypeCaller
8) Somatic variant detection using GATK Mutect2

Working directory structure for running NGS analysis:

1) Paired-end reads (ending in _R1.fastq and _R1.fastq, say A_R1.fastq and A_R2.fastq) in ‘FastQ’ directory
2) CallVars.yml
3) Snakefile
4) 1000G_phase1.indels.hg19.sites.vcf
5) 1000G_phase1.indels.hg19.sites.vcf.idx
6) dbSNP_20180423.vcf
7) dbSNP_20180423.vcf.idx
8) GNOMAD_hg19.vcf
9) GNOMAD_hg19.vcf.idx
10) Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
11) Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
12) Reference Genome Directory ‘UCSCWholeGenomeFasta’ containing below files

  	•genome.dict, genome.fa, genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.fai, genome.fa.pac, genome.fa.sa, GenomeSize.xml

Steps (on the Linux command line interface) for running NGS analysis.

1)	Check the working directory and FastQ files: 

Make sure you are in the working directory that contains all the needed files and folders for running snakemake.

Make sure your Fastq files are in ‘FastQ’ directory and they end in ‘_R1.fastq’ and ‘_R1.fastq’, say A_R1.fastq and A_R2.fastq. 

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



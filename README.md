# CallVars: 

CallVars is an automated, reproducible Snakemake workflow which takes paired-end FastQ files directly to a filtered list of high confidence variants for clinical review. This workflow largely follows [Broad Institute's Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) guidelines for germline short variant discovery (SNPs + Indels) for single sample and also reports a filtered list of somatic variants. 

CallVars can be helpful to anyone working with targeted cancer/rare disease gene panels to find and report variants of clinical relevance. If you think CallVars can help with your study feel free to DM me on twitter (@IAnalyzeGenomes). Feedback/comments/bug reports/contributions are welcome for improvement of this workflow.

# CallVars workflow overview:
CallVars sequentially performs below steps of Next-Gen Sequencing (NGS) analysis.
1) Pre-processing using Cutadapt

	Pre-processing prepares the data for NGS analysis. When DNA or RNA molecules are sequenced using Illumina short reads technology, the machine may sequence into the adapter ligated to the 3’ end of each molecule during library preparation. Consequently, the reads that are output contain the sequence of the molecule of interest and also the adapter sequence. Also, with Illumina sequencing machines, the quality of reads is high at the beginning but degrades towards the 3’ end of the reads. 
	
	CallVars uses Cutadapt to remove adapters from sequencing reads. Cutadapt also trims the read ends with quality below 20 and removes the ambiguous bases (N’s) from the reads ends. 

	You will need to edit the snakefile (line 12 and line 25) to remove below adapters and include adapters used by your library preparation kit.
	
	AGATCGGAAGAGC 
	
	AGATCGGAAGAGC
	 

2) Mapping using BWA
	
	Once the high quality reads are obtained from pre-processing, the next step is mapping them to human reference genome. CallVars use industry standard BWA-mem to map short Illumina paired-end reads to hg19 version of human reference genome. This step generate a Binary Alignment Map also called a BAM file. The reference genome files needed for the analysis were downloaded in 2bit format using below link.
	
	http://hgdownload.cse.ucsc.edu/gbdb/hg19/

3) Sorting using samtools
	
	Now that we have a BAM file, we need to index it. All BAM files need an index, as they tend to be large and the index allows us to perform computationally complex operations on these files without it taking days to complete. Before we index the BAM file we need to sort them by position and remove duplicates. This step performs sorting the BAM file by position.
	
4) Removing duplicates using GATK MarkDuplicates
	
	CallVars uses this tool to locate and remove duplicate reads in a BAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. The MarkDuplicates tool works by comparing sequences in the 5 prime positions of read-pairs in a BAM file.
	
5) Indexing using samtools
	
	CallVars now perform indexing discussed in step 3 using samtools. 
	
6) Base quality score recalibration using GATK BaseRecalibrator and ApplyBQSR
	
	CallVars uses GATK BaseRecalibrator and ApplyBQSR to perform Base Quality Score Recalibration (BQSR). Every base sequenced by machine is assigned a base quality score. These scores play an important role in variant detection. Unfortunately, the scores produced by the machines are subject to various sources of systematic (non-random) technical error. Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment.

	BQSR uses a machine learning approach to model and correct systematic base scoring errors in particular regions of the genome such as homopolymer runs. It works in a two-pass manner, first building a model over all bases in the dataset as well as a set of known variants and writing the model to a recalibration table, as performed by GATK BaseRecalibrator. The second pass actually applies the learned model to correct per-base alignment quality scores to output a recalibrated BAM, as performed by GATK AppyBQSR. 
	

7) Germline variant detection using GATK HaplotypeCaller

	CallVars uses GATK HaplotypeCaller to call germline SNPs and indels via local re-assembly of haplotypes. The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.

8) Somatic variant detection using GATK Mutect2
	
	CallVars uses GATK Mutect2 to call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations. The caller uses a Bayesian somatic genotyping model and uses the assembly-based machinery of HaplotypeCaller. 

9) Functional annotation for germline variants using GATK Funcotator
	
	CallVars uses GATK Funcotator (FUNCtional annOTATOR) to analyze given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file. This tool is a functional annotation tool that allows a user to add annotations to called variants based on a set of data sources, each with its own matching criteria.
	
	Below listed data souces were used for annotation of variants.
	
	Genecode
	Clinvar
	Gnomad
	
	Data source was downloaded using below link.
	https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator --> funcotator_dataSources.v1.6.20190124s.tar.gz
	
10) Functional annotation for somatic variants using GATK Funcotator
	
	This step performs functional annotation as discussed in step 10 for somatic variants.

11) Variant filtration for germline variants using GATK VariantFiltration
	
	 [GATK guidelines](https://software.broadinstitute.org/gatk/documentation/article.php?id=6925) were used to apply generic hard-filtering to add PASS/FAIL tags to variants. Variants are not filtered out based on their PASS/FAIL status.
	 CallVars currently uses gnomAD allele frequency as a key filter to filter and report variants having either genomes or exomes allele frequency less than 0.5% for clinical review. 
	
12) Variant filtration for somatic variants using GATK VariantFiltration 
	
	This step performs filtering as discussed in step 11 for Somatic variants.


# Setting up a working directory to run CallVars:
All the below listed files/folders must be present in the working directory before you run CallVars. Make sure the names of files/folders match exctly as listed.  

	- "FastQ" folder containing paired-end reads ending in _R1.fastq and _R2.fastq (test files A_R1.fastq and A_R2.fastq attached in repo)
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
	The utility program, twoBitToFa (available from the kent src tree), can be used to extract .fa file(s) from this file.  A pre-compiled version of the command line tool can be found at: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
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
		Ensure you run the below command's in working directory.

		Dry run:
		snakemake -np CallVars/Reports/{your_sample_name}_Germline.txt CallVars/Reports/{your_sample_name}_Somatic.txt --cores N
		Real run:
		snakemake CallVars/Reports/{your_sample_name}_Germline.txt CallVars/Reports/{your_sample_name}_Somatic.txt --cores N

For instance, if the names if the FastQ files are A_R1.fastq and A_R2.fastq and number of cores available are 6 then run below command.

		Dry run:
		snakemake -np CallVars/Reports/A_Germline.txt CallVars/Reports/A_Somatic.txt --cores 6
		Real run:
		snakemake CallVars/Reports/A_Germline.txt CallVars/Reports/A_Somatic.txt --cores 6

After the worklow has run successfully, below listed files will be available for clinical review.

	1] CallVars/Reports/A_Germline.txt containing a filtered list of germline variants.
	2] CallVars/Reports/A_Somatic.txt containing a filtered list of somatic variants.
	3] CallVars/Reports/A_Germline_All.vcf containing a full list of germline variants.
	4] CallVars/Reports/A_Somatic_All.vcf containing a full list of somatic variants.


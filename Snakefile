######################################################################################################################################
############################################BEFORE YOU RUN CALLVARS##################################################################
######################################################################################################################################
#1. Make sure to correctly organize the working directory. Read "Setting up a working directory to run CallVars:" section in README.md
#2. Make sure to check the config.yaml file for desired samples or parameters to be used for analysis.
######################################################################################################################################	

configfile: "config.yaml"

#Rule All
rule all:
	input:
		FINAL1=expand("CallVars/Reports/{sample}_Somatic.txt", sample=config["SAMPLE"]),
		FINAL2=expand("CallVars/Reports/{sample}_Germline.txt", sample=config["SAMPLE"])

#Rule to perform adapter trimming using CUTADAPT for files ending in .fastq.gz

rule CUTADAPT_Trim1:
	input:
		FWD="FastQ/{sample}_R1.fastq.gz",
		REV="FastQ/{sample}_R2.fastq.gz"
	output:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"	
	log:
		"CallVars/Logs/{sample}_CUTADAPT-Trimming.log"
	params:
		adp1=expand("{adp1}", adp1=config["adapter1_config"]),
		adp2=expand("{adp2}", adp2=config["adapter2_config"]),
		Q1=expand("{Q1}", Q1=config["q1_config"]),
		Q2=expand("{Q2}", Q2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	run:
		cores = params.cutadapt_threads
		if cores == ['1']:
			print("You provided single CPU core for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} &>{log}")
		else:
			print("You provided ", cores, "CPU cores for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}")

#Rule to perform adapter trimming using CUTADAPT for files ending in .fastq

rule CUTADAPT_Trim2:
	input:
		FWD="FastQ/{sample}_R1.fastq",
		REV="FastQ/{sample}_R2.fastq"
	output:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"
	log:
		"CallVars/Logs/{sample}_CUTADAPT-Trimming.log"
	params:
		adp1=expand("{adp1}", adp1=config["adapter1_config"]),
		adp2=expand("{adp2}", adp2=config["adapter2_config"]),
		Q1=expand("{Q1}", Q1=config["q1_config"]),
		Q2=expand("{Q2}", Q2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	run:	
		cores = params.cutadapt_threads
		if cores == ['1']:
			print("You provided single CPU core for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} &>{log}")
		else:
			print("You provided ", cores, "CPU cores for cutadapt to run. ")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}")

#Rule to mapping of reads to reference genome

rule BWA_Mapping:
	input:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"
	output:
		"CallVars/MappedReads/{sample}.bam"
	params:
		rg=r"@RG\tID:{sample}\tLB:TS1E\tPL:Illumina\tPU=NextSeq550\tSM:{sample}",
		bwa_threads=expand("{bwa_threads}", bwa_threads=config["bwa_threads_config"]),
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"]),
		REF=expand("{REF}", REF=config["Reference"])
	log:
		"CallVars/Logs/{sample}_BWA-Mapping.log"
	shell:
		"bwa mem -R '{params.rg}' -t {params.bwa_threads}  {params.REF} {input.FWD_TRIM} {input.REV_TRIM} | samtools view -Sb - > {output} -@ {params.samtools_threads} 2>{log}"

#Rule to sort the BAM file by position using SAMTOOLS

rule SAMTOOLS_Sort:
	input:
		"CallVars/MappedReads/{sample}.bam"
	output:
		"CallVars/SortedReads/{sample}.bam"
	params:
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"])
	log:
		"CallVars/Logs/{sample}_SAMTOOLS-sort.log"
	shell:
		"samtools sort -T CallVars/SortedReads/{wildcards.sample} "
		"-O bam {input} > {output} -@ {params.samtools_threads} 2>{log}"

#Rule to remove duplicate reads using GATK MarkDuplicates

rule GATK_MarkDuplicates:
	input:
		"CallVars/SortedReads/{sample}.bam"
	output:
		first="CallVars/NoDupReads/{sample}.bam",
		second="/home/amit/Desktop/Shared/BWA-GATK/CallVars/NoDupReads/{sample}_marked_dup_metrics.txt"
	log:
		"CallVars/Logs/{sample}_GATK-MarkDuplicates.log"
	shell:
		"gatk MarkDuplicates --REMOVE_DUPLICATES -I={input} -O={output.first} -M={output.second} --VALIDATION_STRINGENCY=SILENT 2>{log}"

#Rule to index the BAM file using SAMTOOLS

rule SAMTOOLS_Index:
	input:
		"CallVars/NoDupReads/{sample}.bam"
	output:
		"CallVars/NoDupReads/{sample}.bam.bai"
	params:
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"])
	shell:
		"samtools index {input} -@ {params.samtools_threads}"

#Rule to perform part one base quality score recalibration using GATK BaseRecalibrator

rule GATK_BaseRecalibrator:
	input:
		"CallVars/NoDupReads/{sample}.bam"
	output:
		"CallVars/BQSR/{sample}_recal_data.table"
	params:
		 mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		 REF=expand("{REF}", REF=config["Reference"]),
		 dbSNP=expand("{dbSNP}", dbSNP=config["dbSNP_Database"]),
		 MILLS=expand("{MILLS}", MILLS=config["MILLS_Database"]),
		 dbSNPIndels=expand("{dbSNPIndels}", dbSNPIndels=config["dbSNPIndels_Database"])
	log:
		"CallVars/Logs/{sample}_GATK-BaseRecalibrator.log"
	shell:
		"gatk --java-options '{params.mem}' BaseRecalibrator -I {input} -R {params.REF} --known-sites {params.dbSNP} --known-sites {params.MILLS} --known-sites {params.dbSNPIndels} -O {output} &>{log}"

#Rule to perform part two of base quality score recalibration using GATK BQSR

rule GATK_BQSR:
	input:
		BAM="CallVars/NoDupReads/{sample}.bam",
		RECAL="CallVars/BQSR/{sample}_recal_data.table"
	output:
		"CallVars/BQSR/{sample}.bam"
	log:
		"CallVars/Logs/{sample}_GATK-BQSR.log"
	params:
		 mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		 REF=expand("{REF}", REF=config["Reference"]),	
	shell:
		"gatk --java-options '{params.mem}' ApplyBQSR -R {params.REF} -I {input.BAM} --bqsr-recal-file {input.RECAL} -O {output} &>{log}"

#Rule to call germline variants using GATK HaplotypeCaller

rule GATK_HaplotypeCaller:
	input:
		"CallVars/BQSR/{sample}.bam",
	output:
		"CallVars/VCF/{sample}_germline.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-HaplotypeCaller.log"
	params:
		mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		REF=expand("{REF}", REF=config["Reference"]),
		dbSNP=expand("{dbSNP}", dbSNP=config["dbSNP_Database"]),
		TARGET=expand("{TARGET}", TARGET=config["TARGET_config"])
	shell:
		#"gatk --java-options '{params.mem}' HaplotypeCaller -R {input.REF} -I {input.BAM} --dbsnp {input.SNP} -O {output} &>{log}"
		"gatk --java-options '{params.mem}' HaplotypeCaller -R {params.REF} -I {input} --dbsnp {params.dbSNP} -L {params.TARGET} -O {output} &>{log}"

#Rule to call somatic variants using GATK Mutect2
	
rule GATK_Mutect2:
	input:
		"CallVars/BQSR/{sample}.bam",
	output:
		"CallVars/VCF/{sample}_somatic.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Mutect2.log"
	params:
		mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		TARGET=expand("{TARGET}", TARGET=config["TARGET_config"]),
		REF=expand("{REF}", REF=config["Reference"]),
		gnomAD=expand("{REF}", REF=config["GNOMAD_Database"])
	shell:
		#shell("gatk --java-options '{params.mem}' Mutect2 -R {input.REF} -I {input.BAM} -tumor {wildcards.sample} --germline-resource {input.GNOMAD} -O {output} &>{log}")
		"""
		gatk --java-options '{params.mem}' Mutect2 -R {params.REF} -I {input} -L {params.TARGET} -tumor {wildcards.sample} --germline-resource {params.gnomAD} -O {output} &>{log}
		rm CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R1.fastq.gz CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R2.fastq.gz
		rm CallVars/MappedReads/{wildcards.sample}.bam CallVars/SortedReads/{wildcards.sample}.bam CallVars/NoDupReads/{wildcards.sample}.bam
		"""

#Rule to functionally annotate germline variants using GATK Funcotator 

rule GATK_Funcotator_Germline:
	input:
		"CallVars/VCF/{sample}_germline.vcf"
	output:
		"CallVars/VCF/{sample}_germline_func.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Funcotator_Germline.log"
	params:
		 mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		 REF=expand("{REF}", REF=config["Reference"])
	shell:
		"gatk --java-options '{params.mem}' Funcotator -R {params.REF} -V {input} -O {output} --output-file-format VCF --data-sources-path dataSourcesFolder/ --ref-version hg19 &>{log}"

#Rule to functionally annotate somatic variants using GATK Funcotator 

rule GATK_Funcotator_Somatic:
	input:
		"CallVars/VCF/{sample}_somatic.vcf"
	output:
		"CallVars/VCF/{sample}_somatic_func.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Funcotator_Somatic.log"
	params:
		 mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		 REF=expand("{REF}", REF=config["Reference"])
	shell:
		"gatk --java-options '{params.mem}' Funcotator -R {params.REF} -V {input} -O {output} --output-file-format VCF --data-sources-path dataSourcesFolder/ --ref-version hg19 &>{log}"

#Rule to add PASS/FAIL tags to germline variants using GATK VariantFiltration and then filter variants with gnomAD genomes or exomes allele freq less than 0.5%

rule GATK_VariantFiltration_Germline:
	input:
		"CallVars/VCF/{sample}_germline_func.vcf"
	output:
		FILTER="CallVars/Reports/{sample}_Germline_All.vcf",
		ONE="CallVars/TempFiles/{sample}_CallVars_Germline_TempOne.txt",
		TWO="CallVars/TempFiles/{sample}_CallVars_Germline_TempTwo.txt",
		THREE="CallVars/TempFiles/{sample}_CallVars_Germline_TempThree.txt",
		FINALTEMP="CallVars/TempFiles/{sample}_CallVars_Germline_FinalTemp.txt",
		FINAL="CallVars/Reports/{sample}_Germline.txt"
	log:
		"CallVars/Logs/{sample}_GATK-VariantFiltration_Germline.log"
	params:
		mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		REF=expand("{REF}", REF=config["Reference"]),
		gnomAD_Filter=expand("{gnomAD_Filter}", gnomAD_Filter=config["gnomAD_Filter_config"]),
		QD_Filter=expand("{QD_Filter}", QD_Filter=config["QD_Filter_config"]),
		FS_Filter=expand("{FS_Filter}", FS_Filter=config["FS_Filter_config"]),
		MQ_Filter=expand("{MQ_Filter}", MQ_Filter=config["MQ_Filter_config"]),
		MQRankSum_Filter=expand("{MQRankSum_Filter}", MQRankSum_Filter=config["MQRankSum_Filter_config"]),
		ReadPosRankSum_Filter=expand("{ReadPosRankSum_Filter}", ReadPosRankSum_Filter=config["ReadPosRankSum_Filter_config"]),
		SOR_Filter=expand("{SOR_Filter}", SOR_Filter=config["SOR_Filter_config"])
	shell:
		"""
		gatk --java-options '{params.mem}' VariantFiltration -R {params.REF} -O {output.FILTER} -V {input} --filter-expression \"(QD < {params.QD_Filter}) || (FS > {params.FS_Filter}) || (MQ < {params.MQ_Filter}) || (MQRankSum < {params.MQRankSum_Filter}) || (ReadPosRankSum < {params.ReadPosRankSum_Filter}) || (SOR > {params.SOR_Filter})\" --filter-name \"Fail\" &>{log}
		cat {output.FILTER}| cut -f1-7 | tail -n +60 > {output.ONE} || true
		cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION\=\["  '{{ print $2 }}'  | awk -F"|" '{{ print $1"\t"$6"\t"$14"\t"$17"\t"$19"\t"$27"\t"$28*100"\t"$68*100}}' | tail -n +60 > {output.TWO} ||true
		cat {output.FILTER}| cut -f9,10 | tail -n +60 > {output.THREE} || true
		paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
		awk -F"\\t" '{{ if ($14 <= {params.gnomAD_Filter} || $15 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
		"""

#Rule to add PASS/FAIL tags to somatic variants using GATK VariantFiltration and then filter variants with gnomAD genomes or exomes allele freq less than 0.5%

rule GATK_VariantFiltration_Somatic:
	input:
		"CallVars/VCF/{sample}_somatic_func.vcf"
	output:
		FILTER="CallVars/Reports/{sample}_Somatic_All.vcf",
		ONE="CallVars/TempFiles/{sample}_CallVars_Somatic_TempOne.txt",
		TWO="CallVars/TempFiles/{sample}_CallVars_Somatic_TempTwo.txt",
		THREE="CallVars/TempFiles/{sample}_CallVars_Somatic_TempThree.txt",
		FINALTEMP="CallVars/TempFiles/{sample}_CallVars_Somatic_FinalTemp.txt",
		FINAL="CallVars/Reports/{sample}_Somatic.txt"
	log:
		"CallVars/Logs/{sample}_GATK-Funcotator_Somatic.log"
	params:
		mem=expand("{mem}", mem=config["GATK_JAVA_config"]),
		REF=expand("{REF}", REF=config["Reference"]),
		gnomAD_Filter=expand("{gnomAD_Filter}", gnomAD_Filter=config["gnomAD_Filter_config"]),
		QD_Filter=expand("{QD_Filter}", QD_Filter=config["QD_Filter_config"]),
		FS_Filter=expand("{FS_Filter}", FS_Filter=config["FS_Filter_config"]),
		MQ_Filter=expand("{MQ_Filter}", MQ_Filter=config["MQ_Filter_config"]),
		MQRankSum_Filter=expand("{MQRankSum_Filter}", MQRankSum_Filter=config["MQRankSum_Filter_config"]),
		ReadPosRankSum_Filter=expand("{ReadPosRankSum_Filter}", ReadPosRankSum_Filter=config["ReadPosRankSum_Filter_config"]),
		SOR_Filter=expand("{SOR_Filter}", SOR_Filter=config["SOR_Filter_config"])
	shell:
		"""
		gatk --java-options '{params.mem}' VariantFiltration -R {params.REF} -O {output.FILTER} -V {input} --filter-expression \"(QD < {params.QD_Filter}) || (FS > {params.FS_Filter}) || (MQ < {params.MQ_Filter}) || (MQRankSum < {params.MQRankSum_Filter}) || (ReadPosRankSum < {params.ReadPosRankSum_Filter}) || (SOR > {params.SOR_Filter})\" --filter-name \"Fail\" &>{log}
		cat {output.FILTER}| cut -f1-7 | tail -n +71 > {output.ONE} || true
		cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION\=\["  '{{ print $2 }}'  | awk -F"|" '{{ print $1"\t"$6"\t"$14"\t"$17"\t"$19"\t"$27"\t"$28*100"\t"$68*100}}' | tail -n +71 > {output.TWO} ||true
		cat {output.FILTER}| cut -f9,10 | tail -n +71 > {output.THREE} || true
		paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
		awk -F"\\t" '{{ if ($14 <= {params.gnomAD_Filter} || $15 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
		"""
onsuccess:
	print("SUCCESS: CallVars has been run succesfully!!")
onerror:
	print("ERROR!! CallVars has been aborted!! Before you run CallVars make sure to")
	print ("1. correctly organize the working directory. Read \"Setting up a working directory to run CallVars:\" section in README.md")
	print ("2. check the config.yaml file for desired samples or parameters to be used for analysis.")

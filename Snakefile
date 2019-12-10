configfile: "config.yaml"
rule all:
	input:
		FINAL1="CallVars/Reports/A_Somatic.txt",
		FINAL2="CallVars/Reports/A_Germline.txt"

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
		Q1=expand("{adp2}", adp2=config["q1_config"]),
		Q2=expand("{adp2}", adp2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	shell:
		"cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}"

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
		Q1=expand("{adp2}", adp2=config["q1_config"]),
		Q2=expand("{adp2}", adp2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	shell:
		"cutadapt -q {params.Q2},{params.Q1} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}"

rule BWA_Mapping:
	input:
		REF="UCSCWholeGenomeFasta/genome.fa",
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"
	output:
		"CallVars/MappedReads/{sample}.bam"
	params:
		rg=r"@RG\tID:{sample}\tLB:TS1E\tPL:Illumina\tPU=NextSeq550\tSM:{sample}",
		bwa_threads=expand("{bwa_threads}", bwa_threads=config["bwa_threads_config"]),
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"])
	log:
		"CallVars/Logs/{sample}_BWA-Mapping.log"
	shell:
		"bwa mem -R '{params.rg}' -t {params.bwa_threads}  {input.REF} {input.FWD_TRIM} {input.REV_TRIM} | samtools view -Sb - > {output} -@ {params.samtools_threads} 2>{log}"

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

rule SAMTOOLS_Index:
	input:
		"CallVars/NoDupReads/{sample}.bam"
	output:
		"CallVars/NoDupReads/{sample}.bam.bai"
	params:
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"])
	shell:
		"samtools index {input} -@ {params.samtools_threads}"

rule GATK_BaseRecalibrator:
	input:
		BAM="CallVars/NoDupReads/{sample}.bam",
		REF="UCSCWholeGenomeFasta/genome.fa",
		SNP="dbSNP_20180423.vcf",
		MILLS="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
		GNOM="1000G_phase1.indels.hg19.sites.vcf"
	output:
		"CallVars/BQSR/{sample}_recal_data.table"
	params:
		 mem="-Xmx30g -Xmx20g"
	log:
		"CallVars/Logs/{sample}_GATK-BaseRecalibrator.log"
	shell:
		"gatk --java-options '{params.mem}' BaseRecalibrator -I {input.BAM} -R {input.REF} --known-sites {input.SNP} --known-sites {input.MILLS} --known-sites {input.GNOM} -O {output} &>{log}"

rule GATK_BQSR:
	input:
		BAM="CallVars/NoDupReads/{sample}.bam",
		REF="UCSCWholeGenomeFasta/genome.fa",
		RECAL="CallVars/BQSR/{sample}_recal_data.table"
	output:
		"CallVars/BQSR/{sample}.bam"
	log:
		"CallVars/Logs/{sample}_GATK-BQSR.log"
	params:
		 mem="-Xmx30g -Xmx20g"	
	shell:
		"gatk --java-options '{params.mem}' ApplyBQSR -R {input.REF} -I {input.BAM} --bqsr-recal-file {input.RECAL} -O {output} &>{log}"

rule GATK_HaplotypeCaller:
	input:
		BAM="CallVars/BQSR/{sample}.bam",
		REF="UCSCWholeGenomeFasta/genome.fa",
		SNP="dbSNP_20180423.vcf"
	output:
		"CallVars/VCF/{sample}_germline.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-HaplotypeCaller.log"
	params:
		mem="-Xmx30g -Xmx20g",
		TARGET=expand("{TARGET}", TARGET=config["TARGET_config"])
	shell:
		#"gatk --java-options '{params.mem}' HaplotypeCaller -R {input.REF} -I {input.BAM} --dbsnp {input.SNP} -O {output} &>{log}"
		"gatk --java-options '{params.mem}' HaplotypeCaller -R {input.REF} -I {input.BAM} --dbsnp {input.SNP} -L {params.TARGET} -O {output} &>{log}"
		
rule GATK_Mutect2:
	input:
		BAM="CallVars/BQSR/{sample}.bam",
		REF="UCSCWholeGenomeFasta/genome.fa",
		GNOMAD="GNOMAD_hg19.vcf",
		TARGET="TARGET.bed"
	output:
		"CallVars/VCF/{sample}_somatic.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Mutect2.log"
	params:
		mem="-Xmx30g -Xmx20g",
		TARGET=expand("{TARGET}", TARGET=config["TARGET_config"])
	run:
		#shell("gatk --java-options '{params.mem}' Mutect2 -R {input.REF} -I {input.BAM} -tumor {wildcards.sample} --germline-resource {input.GNOMAD} -O {output} &>{log}")
		shell("gatk --java-options '{params.mem}' Mutect2 -R {input.REF} -I {input.BAM} -L {input.TARGET} -tumor {wildcards.sample} --germline-resource {input.GNOMAD} -O {output} &>{log}")	
		shell("rm CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R1.fastq.gz CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R2.fastq.gz")
		shell("rm CallVars/MappedReads/{wildcards.sample}.bam CallVars/SortedReads/{wildcards.sample}.bam CallVars/NoDupReads/{wildcards.sample}.bam")

rule GATK_Funcotator_Germline:
	input:
		REF="UCSCWholeGenomeFasta/genome.fa",
		VCF="CallVars/VCF/{sample}_germline.vcf"
	output:
		"CallVars/VCF/{sample}_germline_func.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Funcotator_Germline.log"
	params:
		 mem="-Xmx30g -Xmx20g"
	run:
		shell("gatk --java-options '{params.mem}' Funcotator -R {input.REF} -V {input.VCF} -O {output} --output-file-format VCF --data-sources-path dataSourcesFolder/ --ref-version hg19 &>{log}")

rule GATK_Funcotator_Somatic:
	input:
		REF="UCSCWholeGenomeFasta/genome.fa",
		VCF="CallVars/VCF/{sample}_somatic.vcf"
	output:
		"CallVars/VCF/{sample}_somatic_func.vcf"
	log:
		"CallVars/Logs/{sample}_GATK-Funcotator_Somatic.log"
	params:
		 mem="-Xmx30g -Xmx20g"
	run:
		shell("gatk --java-options '{params.mem}' Funcotator -R {input.REF} -V {input.VCF} -O {output} --output-file-format VCF --data-sources-path dataSourcesFolder/ --ref-version hg19 &>{log}")


rule GATK_VariantFiltration_Germline:
	input:
		REF="UCSCWholeGenomeFasta/genome.fa",
		VCF="CallVars/VCF/{sample}_germline_func.vcf"
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
		mem="-Xmx30g -Xmx20g",
		gnomAD_Filter=expand("{gnomAD_Filter}", gnomAD_Filter=config["gnomAD_Filter_config"]),
		QD_Filter=expand("{QD_Filter}", QD_Filter=config["QD_Filter_config"]),
		FS_Filter=expand("{FS_Filter}", FS_Filter=config["FS_Filter_config"]),
		MQ_Filter=expand("{MQ_Filter}", MQ_Filter=config["MQ_Filter_config"]),
		MQRankSum_Filter=expand("{MQRankSum_Filter}", MQRankSum_Filter=config["MQRankSum_Filter_config"]),
		ReadPosRankSum_Filter=expand("{ReadPosRankSum_Filter}", ReadPosRankSum_Filter=config["ReadPosRankSum_Filter_config"]),
		SOR_Filter=expand("{SOR_Filter}", SOR_Filter=config["SOR_Filter_config"])
	shell:
		"""
		gatk --java-options '{params.mem}' VariantFiltration -R {input.REF} -O {output.FILTER} -V {input.VCF} --filter-expression \"(QD < {params.QD_Filter}) || (FS > {params.FS_Filter}) || (MQ < {params.MQ_Filter}) || (MQRankSum < {params.MQRankSum_Filter}) || (ReadPosRankSum < {params.ReadPosRankSum_Filter}) || (SOR > {params.SOR_Filter})\" --filter-name \"Fail\" &>{log}
		cat {output.FILTER}| cut -f1-7 | tail -n +60 > {output.ONE} || true
		cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION\=\["  '{{ print $2 }}'  | awk -F"|" '{{ print $1"\t"$6"\t"$14"\t"$17"\t"$19"\t"$27"\t"$28*100"\t"$68*100}}' | tail -n +60 > {output.TWO} ||true
		cat {output.FILTER}| cut -f9,10 | tail -n +60 > {output.THREE} || true
		paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
		awk -F"\\t" '{{ if ($14 <= {params.gnomAD_Filter} || $15 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
		"""

rule GATK_VariantFiltration_Somatic:
	input:
		REF="UCSCWholeGenomeFasta/genome.fa",
		VCF="CallVars/VCF/{sample}_somatic_func.vcf"
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
		mem="-Xmx30g -Xmx20g",
		gnomAD_Filter=expand("{gnomAD_Filter}", gnomAD_Filter=config["gnomAD_Filter_config"]),
	shell:
		"""
		gatk --java-options '{params.mem}' VariantFiltration -R {input.REF} -O {output.FILTER} -V {input.VCF} --filter-expression \"(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 3.0)\" --filter-name \"Fail\" &>{log}
		cat {output.FILTER}| cut -f1-7 | tail -n +71 > {output.ONE} || true
		cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION\=\["  '{{ print $2 }}'  | awk -F"|" '{{ print $1"\t"$6"\t"$14"\t"$17"\t"$19"\t"$27"\t"$28*100"\t"$68*100}}' | tail -n +71 > {output.TWO} ||true
		cat {output.FILTER}| cut -f9,10 | tail -n +71 > {output.THREE} || true
		paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
		awk -F"\\t" '{{ if ($14 <= {params.gnomAD_Filter} || $15 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
		"""
onsuccess:
	print("SUCCESS: CallVars has been run succesfully!!")
onerror:
	print("ERROR!! CallVars has been aborted!!")


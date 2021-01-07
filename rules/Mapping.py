#Rule to mapping of reads to reference genome

rule BWA_Mapping:
	input:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"
	output:
		"CallVars/MappedReads/{sample}.bam"
	params:
		rg=r"-R '@RG\tID:{sample}\tSM:{sample}'",
		bwa_threads=expand("{bwa_threads}", bwa_threads=config["bwa_threads_config"]),
		samtools_threads=expand("{samtools_threads}", samtools_threads=config["samtools_threads_config"]),
		REF=expand("{REF}", REF=config["Reference"])
	log:
		"CallVars/Logs/{sample}_BWA-Mapping.log"
	shell:
		"bwa mem {params.rg} -t {params.bwa_threads}  {params.REF} {input.FWD_TRIM} {input.REV_TRIM} | samtools view -S -b > {output} -@ {params.samtools_threads} 2>{log}"

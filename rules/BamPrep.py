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
		second="CallVars/NoDupReads/{sample}_marked_dup_metrics.txt"
	log:
		"CallVars/Logs/{sample}_GATK-MarkDuplicates.log"
	shell:
		gatk MarkDuplicates --REMOVE_DUPLICATES true  -I {input} -O {output.first} -M {output.second}  2>{log}
     
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

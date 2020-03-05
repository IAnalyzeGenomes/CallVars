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
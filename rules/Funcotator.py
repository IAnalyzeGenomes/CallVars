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
		 REF=expand("{REF}", REF=config["Reference"]),
		 DataSource=expand("{DataSource}", DataSource=config["DataSource_config"]),
		 RefSource=expand("{RefSource}", RefSource=config["Ref_config"])
	shell:
		"gatk --java-options '{params.mem}' Funcotator -R {params.REF} -V {input} -O {output} --output-file-format VCF --data-sources-path '{params.DataSource}'/ --ref-version '{params.RefSource}' &>{log}"

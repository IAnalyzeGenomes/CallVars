# Rule to call germline variants using GATK HaplotypeCaller

rule GATK_HaplotypeCaller:
    input:
        BAM = "CallVars/BQSR/{sample}.bam",
        # COV1 = "CallVars/BQSR/{sample}_PerBaseCov.txt",
        # COV2 = "CallVars/BQSR/{sample}_PerBaseCov_LessThan20.txt"
    output:
        "CallVars/VCF/{sample}_germline.vcf"
    log:
        "CallVars/Logs/{sample}_GATK-HaplotypeCaller.log"
    params:
        mem = expand("{mem}", mem=config["GATK_JAVA_config"]),
        REF = expand("{REF}", REF=config["Reference"]),
        dbSNP = expand("{dbSNP}", dbSNP=config["dbSNP_Database"]),
        TARGET = expand("{TARGET}", TARGET=config["TARGET_config"])
    shell:
        "gatk --java-options '{params.mem}' HaplotypeCaller -R {params.REF} -I {input.BAM} --dbsnp {params.dbSNP} -L {params.TARGET} -O {output} &>{log}"

# Rule to call somatic variants using GATK Mutect2

rule GATK_Mutect2:
    input:
        "CallVars/BQSR/{sample}.bam",
    output:
        "CallVars/VCF/{sample}_somatic.vcf"
    log:
        "CallVars/Logs/{sample}_GATK-Mutect2.log"
    params:
        mem = expand("{mem}", mem=config["GATK_JAVA_config"]),
        TARGET = expand("{TARGET}", TARGET=config["TARGET_config"]),
        REF = expand("{REF}", REF=config["Reference"]),
        gnomAD = expand("{REF}", REF=config["GNOMAD_Database"])
    shell:
        """
        gatk --java-options '{params.mem}' Mutect2 -R {params.REF} -I {input} -L {params.TARGET} -tumor {wildcards.sample} --germline-resource {params.gnomAD} -O {output} &>{log}
        #rm CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R1.fastq.gz CallVars/TrimmedReads/{wildcards.sample}_Trimmed_R2.fastq.gz
        #rm CallVars/MappedReads/{wildcards.sample}.bam CallVars/SortedReads/{wildcards.sample}.bam CallVars/NoDupReads/{wildcards.sample}.bam
        """

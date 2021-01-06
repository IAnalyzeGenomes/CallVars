# Rule to call germline variants using GATK HaplotypeCaller

rule GATK_HaplotypeCaller:
    input:
        BAM = "CallVars/NoDupReads/{sample}.bam",
        index = "CallVars/NoDupReads/{sample}.bam.bai"
    output:
        VCF= "CallVars/VCF/{sample}_germline.vcf",
        GATK_BAM = "CallVars/NoDupReads/{sample}_GATK.bam"
    log:
        "CallVars/Logs/{sample}_GATK-HaplotypeCaller.log"
    params:
        mem = expand("{mem}", mem=config["GATK_JAVA_config"]),
        REF = expand("{REF}", REF=config["Reference"]),
        dbSNP = expand("{dbSNP}", dbSNP=config["dbSNP_Database"]),
        TARGET = expand("{TARGET}", TARGET=config["TARGET_config"])
    shell:
        "gatk --java-options '{params.mem}' HaplotypeCaller -R {params.REF} -I {input.BAM} --dbsnp {params.dbSNP} -L {params.TARGET} -O {output.VCF} --bamout {output.GATK_BAM} &>{log}"

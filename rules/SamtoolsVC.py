#Rule for vanriant calling using SAMTOOLS

rule SAMTOOLS_BCF:
    input:
        bam = "CallVars/NoDupReads/{sample}.bam",
        index = "CallVars/NoDupReads/{sample}.bam.bai"
    output:
        "CallVars/Reports/{sample}_samtools.bcf"
    params:
        REF = expand("{REF}", REF=config["Reference"]),
        TARGET = expand("{TARGET}", TARGET=config["TARGET_config"]),
        THREADS= expand("{THREADS}", THREADS=config["samtools_threads_config"])
    log:
        "CallVars/Logs/{sample}_SAMTOOLS-VARIANT-CALL.log"
    shell:
        "samtools mpileup -uf {params.REF} -l {params.TARGET} {input.bam}"
        "| bcftools call --threads {params.THREADS} --ploidy GRCh38 -mv -Ob -o {output} &>{log} || true"
        
rule SAMTOOLS_VCF:
    input:
        "CallVars/Reports/{sample}_samtools.bcf"
    output:
        "CallVars/Reports/{sample}_samtools.vcf"
    shell: 
        "bcftools view {input} > {output}"     




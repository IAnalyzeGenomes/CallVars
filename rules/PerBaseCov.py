# Rule to find per base coverage
rule Bedtools_PerBaseCov:
    input:
        BAM = "CallVars/BQSR/{sample}.bam"
    output:
        COV1 = "CallVars/BQSR/{sample}_PerBaseCov.txt",
        COV2 = "CallVars/BQSR/{sample}_PerBaseCov_LessThan20.txt"
    log:
        "CallVars/Logs/{sample}_Bedtools_PerBaseCov.log"
    params:
        TARGET = expand("{TARGET}", TARGET=config["TARGET_config"])
    shell:
        """
        bedtools coverage -d -abam {input.BAM} -b {params.TARGET} &> {output.COV1}
        awk '$6 < 20' {output.COV1} &> {output.COV2}
        """

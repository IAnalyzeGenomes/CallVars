configfile: "config.yaml"

# Target Rule
rule all:
    input:
        Report1 = expand(
            "CallVars/Reports/{sample}_Somatic.txt", sample=config["SAMPLE"]),
        Report2 = expand(
            "CallVars/Reports/{sample}_Germline.txt", sample=config["SAMPLE"]),
        Report3 = expand(
            "CallVars/Reports/{sample}_samtools.vcf", sample=config["SAMPLE"]),
        Report4 = expand(
            "CallVars/BQSR/{sample}_PerBaseCov.txt", sample=config["SAMPLE"]),
        Report5 = expand(
            "CallVars/BQSR/{sample}_PerBaseCov_LessThan20.txt", sample=config["SAMPLE"])

# Modules
include: "rules/AdapterTrim.py"
include: "rules/Mapping.py"
include: "rules/BamPrep.py"
include: "rules/SamtoolsVC.py"
include: "rules/BQSR.py"
include: "rules/PerBaseCov.py"
include: "rules/VariantCalling.py"
include: "rules/Funcotator.py"
include: "rules/VariantFiltration.py"

onsuccess:
    print("----------------------------------------------")
    print("--SUCCESS: CallVars has been run succesfully--")
    print("----------------------------------------------")
onerror:
    print("#############################################################################################################################")
    print("ERROR!! CallVars has been aborted!! Before you run CallVars make sure to")
    print("1. correctly organize the working directory. Read \"Setting up a working directory to run CallVars:\" section in README.md")
    print("2. check the config.yaml file for desired samples or parameters to be used for analysis.")
    print("#############################################################################################################################")

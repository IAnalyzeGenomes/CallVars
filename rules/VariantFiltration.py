# Rule to add PASS/FAIL tags to germline variants using GATK VariantFiltration and then filter variants with gnomAD genomes or exomes allele freq less than 0.5%

rule GATK_VariantFiltration_Germline:
    input:
        "CallVars/VCF/{sample}_germline_func.vcf"
    output:
        FILTER = "CallVars/Reports/{sample}_Germline_All.vcf",
        ONE = "CallVars/TempFiles/{sample}_CallVars_Germline_TempOne.txt",
        TWO = "CallVars/TempFiles/{sample}_CallVars_Germline_TempTwo.txt",
        THREE = "CallVars/TempFiles/{sample}_CallVars_Germline_TempThree.txt",
        FINALTEMP = "CallVars/TempFiles/{sample}_CallVars_Germline_FinalTemp.txt",
        FINAL = "CallVars/Reports/{sample}_Germline.txt"
    log:
        "CallVars/Logs/{sample}_GATK-VariantFiltration_Germline.log"
    params:
        mem = expand("{mem}", mem=config["GATK_JAVA_config"]),
        REF = expand("{REF}", REF=config["Reference"]),
        REF_version = expand(
            "{REF_version}", REF_version=config["Ref_config"]),
        gnomAD_Filter = expand(
            "{gnomAD_Filter}", gnomAD_Filter=config["gnomAD_Filter_config"]),
        QD_Filter = expand(
            "{QD_Filter}", QD_Filter=config["QD_Filter_config"]),
        FS_Filter = expand(
            "{FS_Filter}", FS_Filter=config["FS_Filter_config"]),
        MQ_Filter = expand(
            "{MQ_Filter}", MQ_Filter=config["MQ_Filter_config"]),
        MQRankSum_Filter = expand(
            "{MQRankSum_Filter}", MQRankSum_Filter=config["MQRankSum_Filter_config"]),
        ReadPosRankSum_Filter = expand(
            "{ReadPosRankSum_Filter}", ReadPosRankSum_Filter=config["ReadPosRankSum_Filter_config"]),
        SOR_Filter = expand(
            "{SOR_Filter}", SOR_Filter=config["SOR_Filter_config"]),
        CutGermline = expand(
            "{CutGermline}", CutGermline=config["CutGermline_config"])
    run:
        REF_genome = params.REF_version
        #print(f"The reference genome used is {REF_genome}")
        if (REF_genome == ['hg38']):
            shell("""
					echo "The reference genome is HG38"
					gatk --java-options '{params.mem}' VariantFiltration -R {params.REF} -O {output.FILTER} -V {input} --filter-expression \"(QD < {params.QD_Filter}) || (FS > {params.FS_Filter}) || (MQ < {params.MQ_Filter}) || (MQRankSum < {params.MQRankSum_Filter}) || (ReadPosRankSum < {params.ReadPosRankSum_Filter}) || (SOR > {params.SOR_Filter})\" --filter-name \"Fail\" &>{log}
					cat {output.FILTER}| cut -f1-7 | tail -n +{params.CutGermline} > {output.ONE} || true
					cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION="  '{{ print $2 }}'  | awk -F"|" '{{ print substr($1,2)"\t"$6"\t"$8"\t"$12"\t"$14"\t"$17"\t"$27"\t"$29"\t"$31"\t"$32"\t"$33"\t"$34"\t"$43"\t"$47*100"\t"$87*100}}' | tail -n +{params.CutGermline} > {output.TWO} ||true
					cat {output.FILTER}| cut -f9,10 | tail -n +{params.CutGermline}  > {output.THREE} || true
					paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
					awk -F"\\t" '{{ if ($21 <= {params.gnomAD_Filter} || $22 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
				""")
        else:
            shell("""
					echo "The reference genome is HG19"
					gatk --java-options '{params.mem}' VariantFiltration -R {params.REF} -O {output.FILTER} -V {input} --filter-expression \"(QD < {params.QD_Filter}) || (FS > {params.FS_Filter}) || (MQ < {params.MQ_Filter}) || (MQRankSum < {params.MQRankSum_Filter}) || (ReadPosRankSum < {params.ReadPosRankSum_Filter}) || (SOR > {params.SOR_Filter})\" --filter-name \"Fail\" &>{log}
					cat {output.FILTER}| cut -f1-7 | tail -n +{params.CutGermline} > {output.ONE} || true
					cat {output.FILTER}| cut -f8 | awk -F"FUNCOTATION="  '{{ print $2 }}'  | awk -F"|" '{{ print substr($1,2)"\t"$6"\t"$8"\t"$12"\t"$14"\t"$17"\t"$27"\t"$29"\t"$31"\t"$32"\t"$33"\t"$34"\t"$43"\t"$47*100"\t"$87*100}}' | tail -n +{params.CutGermline} > {output.TWO} ||true
					cat {output.FILTER}| cut -f9,10 | tail -n +{params.CutGermline}  > {output.THREE} || true
					paste {output.ONE} {output.TWO} {output.THREE} > {output.FINALTEMP} || true
					awk -F"\\t" '{{ if ($21 <= {params.gnomAD_Filter} || $22 <= {params.gnomAD_Filter}) {{print}}}}' {output.FINALTEMP} > {output.FINAL} || true
				""")

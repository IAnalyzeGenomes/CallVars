#Rule to perform adapter trimming using CUTADAPT for files ending in .fastq.gz

rule CUTADAPT_Trim1:
	input:
		FWD="FastQ/{sample}_R1.fastq.gz",
		REV="FastQ/{sample}_R2.fastq.gz"
	output:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"	
	log:
		"CallVars/Logs/{sample}_CUTADAPT-Trimming.log"
	params:
		adp1=expand("{adp1}", adp1=config["adapter1_config"]),
		adp2=expand("{adp2}", adp2=config["adapter2_config"]),
		Q1=expand("{Q1}", Q1=config["q1_config"]),
		Q2=expand("{Q2}", Q2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	run:
		cores = params.cutadapt_threads
		if cores == ['1']:
			print("You provided single CPU core for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} &>{log}")
		else:
			print("You provided ", cores, "CPU cores for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}")

#Rule to perform adapter trimming using CUTADAPT for files ending in .fastq

rule CUTADAPT_Trim2:
	input:
		FWD="FastQ/{sample}_R1.fastq",
		REV="FastQ/{sample}_R2.fastq"
	output:
		FWD_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R1.fastq.gz",
		REV_TRIM="CallVars/TrimmedReads/{sample}_Trimmed_R2.fastq.gz"
	log:
		"CallVars/Logs/{sample}_CUTADAPT-Trimming.log"
	params:
		adp1=expand("{adp1}", adp1=config["adapter1_config"]),
		adp2=expand("{adp2}", adp2=config["adapter2_config"]),
		Q1=expand("{Q1}", Q1=config["q1_config"]),
		Q2=expand("{Q2}", Q2=config["q2_config"]),
		cutadapt_threads=expand("{cutadapt_threads}", cutadapt_threads=config["cutadapt_threads_config"])
	run:	
		cores = params.cutadapt_threads
		if cores == ['1']:
			print("You provided single CPU core for cutadapt to run.")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} &>{log}")
		else:
			print("You provided ", cores, "CPU cores for cutadapt to run. ")
			shell("cutadapt -q {params.Q1},{params.Q2} --trim-n -a {params.adp1} -A {params.adp2} -o {output.FWD_TRIM} -p {output.REV_TRIM}  {input} -j {params.cutadapt_threads} &>{log}")
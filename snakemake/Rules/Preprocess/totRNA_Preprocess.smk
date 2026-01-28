########################	
###fastqc on raw data###
########################
rule Fastqc_raw:
	input:
		pair1=f"{Rawdir}/{{raw_sample}}_R1_001.fastq.gz", 
		pair2=f"{Rawdir}/{{raw_sample}}_R2_001.fastq.gz"
	output:
		pair1=f"{Workdir}/FastQC/raw_data/{{raw_sample}}_R1_001_fastqc.zip",
		pair2=f"{Workdir}/FastQC/raw_data/{{raw_sample}}_R2_001_fastqc.zip"
	params:
		outdir=f"{Workdir}/FastQC/raw_data/"
	log:
		f"{Workdir}/logs/FastQC/raw_data/err_fastqc_{{raw_sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate fastqc
		fastqc {input.pair1} -o {params.outdir} 2> {log}
		fastqc {input.pair2} -o {params.outdir} 2>> {log}
		"""

#################		
###Trimmomatic###
#################
rule Trimmomatic:
	input:
		inputfiles
	output:
		temp(f"{Workdir}/trimmed_data/{{sample}}_trimmed_1P.fastq.gz"),
		temp(f"{Workdir}/trimmed_data/{{sample}}_trimmed_2P.fastq.gz"),
		log=f"{Workdir}/logs/trimmomatic/err_trimmomatic_{{sample}}.txt"
	params:
		outdir=f"{Workdir}/trimmed_data/{{sample}}_trimmed.fastq.gz",
		adapters=config["adapters"]
	resources:
		runtime=300,
		mem_mb=10000
	threads: 
		6
	shell:
		"""
		micromamba activate Trimmomatic
		trimmomatic PE -threads {threads} -phred33 {input[0]} {input[1]} -baseout {params.outdir} ILLUMINACLIP:{params.adapters}:2:30:10:8:TRUE LEADING:28 TRAILING:28 SLIDINGWINDOW:4:15 MINLEN:30 2> {output.log}
		"""

############################
###fastqc on trimmed data###
############################
rule Fastqc_trim:
	input:
		pair1 = f"{Workdir}/trimmed_data/{{sample}}_trimmed_1P.fastq.gz",
		pair2 = f"{Workdir}/trimmed_data/{{sample}}_trimmed_2P.fastq.gz"
	output:
		f"{Workdir}/FastQC/trimmed_data/{{sample}}_trimmed_1P_fastqc.zip",
		f"{Workdir}/FastQC/trimmed_data/{{sample}}_trimmed_2P_fastqc.zip"
	params:
		outdir=f"{Workdir}/FastQC/trimmed_data/"
	log:
		f"{Workdir}/logs/FastQC/trimmed_data/err_fastqc_{{sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate fastqc
		fastqc {input.pair1} -o {params.outdir} 2> {log}
		fastqc {input.pair2} -o {params.outdir} 2>> {log}
		"""

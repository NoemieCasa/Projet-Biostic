# ============================================================
# FastQC sur les FASTQ bruts
# ============================================================
rule Fastqc_raw:
	input:
        	r1 = lambda wc: f"{DataPath}/{wc.sample}_R1_001.fastq.gz",
        	r2 = lambda wc: f"{DataPath}/{wc.sample}_R2_001.fastq.gz"
	output:
        	r1_html = f"{Workdir}/results/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
        	r2_html = f"{Workdir}/results/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html"
	params:
		outdir=f"{Workdir}/FastQC/fastqc_raw/"
	log:
		f"{Workdir}/logs/FastQC/results/err_fastqc_{{raw_sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate fastqc
		fastqc {input.r1} -o {params.outdir} 2> {log}
		fastqc {input.r2} -o {params.outdir} 2>> {log}
		"""	

# ============================================================
# Trimming avec Trimmomatic
# ============================================================
rule Trimmomatic:
	input:
        	r1 = f"{DataPath}/{wc.sample}_R1_001.fastq.gz",
        	r2 = f"{DataPath}/{wc.sample}_R2_001.fastq.gz"
	output:
        	r1_paired = f"{Workdir}/results/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        	r1_unpaired = f"{Workdir}/results/Trimming/{{sample}}_R1_001_unpaired.fastq.gz",
        	r2_paired = f"{Workdir}/results/Trimming/{{sample}}_R2_001_paired.fastq.gz",
        	r2_unpaired = f"{Workdir}/results/Trimming/{{sample}}_R2_001_unpaired.fastq.gz"
		log=f"{Workdir}/logs/trimmomatic/err_trimmomatic_{{sample}}.txt"
	params:
		outdir=f"{Workdir}/results/Trimming/{{sample}}_trimmed.fastq.gz",
		adapters=config["adapters"]
	resources:
		runtime=300,
		mem_mb=10000
	threads: 
		6
	shell:
        	
		"""
		micromamba activate trimmomatic
       		mkdir -p {Workdir}/results/Trimming

        	trimmomatic PE \
        	{input.r1} {input.r2} \
        	{output.r1_paired} {output.r1_unpaired} \
        	{output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:{conda_prefix}/share/trimmomatic/adapters/{config[trimmomatic_adapters]}:2:30:10:2:True \
        LEADING:{config[leading]} \
        TRAILING:{config[trailing]} \
        MINLEN:{config[minlen]}
        """
# ============================================================
# FastQC après trimming
# ============================================================
rule Fastqc_trim:
    	input:
        	r1 = f"{Workdir}/results/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        	r2 = f"{Workdir}/results/Trimming/{{sample}}_R2_001_paired.fastq.gz"
    	output:
        	r1_html = f"{Workdir}/results/fastqc_trimmed/{{sample}}/{{sample}}_R1_001_paired_fastqc.html",
        	r2_html = f"{Workdir}/results/fastqc_trimmed/{{sample}}/{{sample}}_R2_001_paired_fastqc.html"
	params:
		outdir=f"{Workdir}/results/fastqc_trimmed/"
	log:
		f"{Workdir}/logs/results/fastqc_trimmed/err_fastqc_{{sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate fastqc
		fastqc {input.pair1} -o {params.outdir} 2> {log}
		fastqc {input.pair2} -o {params.outdir} 2>> {log}
		"""

# ============================================================
# Alignement avec Star
# ============================================================
rule alignement:
	input:


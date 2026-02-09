# ============================================================
# FastQC sur les FASTQ bruts
# ============================================================
rule obtain_fastqc:
    input:
        r1 = lambda wc: f"{config['DataPath']}/{wc.sample}_R1_001.fastq.gz",
        r2 = lambda wc: f"{config['DataPath']}/{wc.sample}_R2_001.fastq.gz"
    output:
        r1_html = f"{config['workdir']}/results/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
        r2_html = f"{config['workdir']}/results/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html"
    threads: config["fastqc_threads"]
    shell:
        """
	micromamba activate fastqc
        mkdir -p {config['workdir']}/results/fastqc_raw/{wildcards.sample}
        fastqc -t {threads} -o {config['workdir']}/results/fastqc_raw/{wildcards.sample} {input.r1}
        fastqc -t {threads} -o {config['workdir']}/results/fastqc_raw/{wildcards.sample} {input.r2}
        """

# ============================================================
# Trimming avec Trimmomatic
# ============================================================
rule trimming:
    input:
        r1 = lambda wc: f"{config['DataPath']}/{wc.sample}_R1_001.fastq.gz",
        r2 = lambda wc: f"{config['DataPath']}/{wc.sample}_R2_001.fastq.gz"
    output:
        r1_paired = f"{config['workdir']}/results/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r1_unpaired = f"{config['workdir']}/results/Trimming/{{sample}}_R1_001_unpaired.fastq.gz",
        r2_paired = f"{config['workdir']}/results/Trimming/{{sample}}_R2_001_paired.fastq.gz",
        r2_unpaired = f"{config['workdir']}/results/Trimming/{{sample}}_R2_001_unpaired.fastq.gz"
    threads: 1
    shell:
        """
	micromamba activate trimmomatic
        mkdir -p {config['workdir']}/results/Trimming

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
rule run_fastqc_trimmed:
    input:
        r1 = f"{config['workdir']}/results/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r2 = f"{config['workdir']}/results/Trimming/{{sample}}_R2_001_paired.fastq.gz"
    output:
        r1_html = f"{config['workdir']}/results/fastqc_trimmed/{{sample}}/{{sample}}_R1_001_paired_fastqc.html",
        r2_html = f"{config['workdir']}/results/fastqc_trimmed/{{sample}}/{{sample}}_R2_001_paired_fastqc.html"
    threads: config["fastqc_threads"]
    shell:
        """
	micromamba activate fastqc
        mkdir -p {config['workdir']}/results/fastqc_trimmed/{wildcards.sample}
        fastqc -t {threads} -o {config['workdir']}/results/fastqc_trimmed/{wildcards.sample} {input.r1}
        fastqc -t {threads} -o {config['workdir']}/results/fastqc_trimmed/{wildcards.sample} {input.r2}
        """

# ============================================================
# Alignement avec Star
# ============================================================
rule alignement:
	input:


# ============================================================
# FastQC sur les FASTQ bruts
# ============================================================
rule Fastqc_raw:
    input:
        r1 = f"{DataPath}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{DataPath}/{{sample}}_R2_001.fastq.gz"
    output:
        r1_html = f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
        r2_html = f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html"
    params:
        outdir=f"{Workdir}/results/fastqc_raw/"
    log:
        f"{Workdir}/logs/fasrqc_raw/err_fastqc_{{sample}}.txt"
    threads: 1
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
        r1 = f"{DataPath}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{DataPath}/{{sample}}_R2_001.fastq.gz"
    output:
        r1_paired = f"{Workdir}/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r1_unpaired = f"{Workdir}/Trimming/{{sample}}_R1_001_unpaired.fastq.gz",
        r2_paired = f"{Workdir}/Trimming/{{sample}}_R2_001_paired.fastq.gz",
        r2_unpaired = f"{Workdir}/Trimming/{{sample}}_R2_001_unpaired.fastq.gz"
    log:
        f"{Workdir}/logs/trimmomatic/err_trimmomatic_{{sample}}.txt"
    params:
        outdir=f"{Workdir}/Trimming/{{sample}}_trimmed.fastq.gz",
        adapters=config["adapters"]
    resources:
        runtime=300,
        mem_mb=10000
    threads: 6
    shell:
        """
        micromamba activate trimmomatic
        mkdir -p {Workdir}/results/Trimming
        trimmomatic PE \
        {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:{params.adapters}:2:30:10:2:True \
        LEADING:{config[leading]} \
        TRAILING:{config[trailing]} \
        MINLEN:{config[minlen]}
        """

# ============================================================
# FastQC après trimming
# ============================================================
rule Fastqc_trim:
    input:
        r1 = f"{Workdir}/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r2 = f"{Workdir}/Trimming/{{sample}}_R2_001_paired.fastq.gz"
    output:
        r1_html = f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R1_001_paired_fastqc.html",
        r2_html = f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R2_001_paired_fastqc.html"
    params:
        outdir=f"{Workdir}/fastqc_trimmed/"
    log:
        f"{Workdir}/logs/fastqc_trimmed/err_fastqc_{{sample}}.txt"
    threads: 1
    shell:
        """
        micromamba activate fastqc
        fastqc {input.r1} -o {params.outdir} 2> {log}
        fastqc {input.r2} -o {params.outdir} 2>> {log}
        """

# ============================================================
# Alignement avec Star
# ============================================================
rule star:
    input:
        index=f"{Workdir}/genome_index/star",
        r1 = f"{Workdir}/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r2 = f"{Workdir}/Trimming/{{sample}}_R2_001_paired.fastq.gz"
    output:
        bam=temp(f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"),
        summary=f"{Workdir}/alignment/star/{{sample}}.starLog.final.out"
    params:
        outpre=f"{Workdir}/alignment/star/{{sample}}.star",
        index=f"{Workdir}/DataBases/HG38_fasta/STAR_index/"
    log:
        f"{Workdir}/logs/star/err_star_{{sample}}.txt"
    resources:
        mem_mb=32000,
        runtime="5h"
    threads: 8
    shell:
        """
        micromamba activate STAR
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.r1}

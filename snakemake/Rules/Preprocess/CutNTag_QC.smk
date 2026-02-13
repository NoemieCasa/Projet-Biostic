# ============================================================
# FastQC sur les FASTQ bruts
# ============================================================
rule Fastqc_raw:
    input:
        unpack(get_raw_fastq)
    output:
        r1_html = f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
        r2_html = f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html"
    params:
        # On crée un dossier par sample pour éviter que FastQC ne mélange tout
        outdir = f"{Workdir}/FastQC/fastqc_raw/{{sample}}"
    log:
        f"{Workdir}/logs/fastqc_raw/err_fastqc_{{sample}}.txt"
    threads: 1
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate fastqc
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} -o {params.outdir} 2> {log}
        """

# ============================================================
# Trimming avec Trimmomatic
# ============================================================
rule Trimmomatic:
    input:
        unpack(get_raw_fastq)
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
        eval "$(micromamba shell hook --shell=bash)"
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
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate fastqc
        fastqc {input.r1} -o {params.outdir} 2> {log}
        fastqc {input.r2} -o {params.outdir} 2>> {log}
        """

# ============================================================
# Alignement avec Star
# ============================================================
rule star:
    input:
        index=f"/home/iguerin2024@ec-nantes.fr/scratch/cutNtag/genome_index/star",
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
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate STAR
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outpre} 2> {log}
        """

# ============================================================
# Samtools stats
# ============================================================
rule samtools_stats:
    input:
        star=f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"
    output:
        stats_star=f"{Workdir}/samstats/star/{{sample}}.star.stats"
    params:
        out_star=f"{Workdir}/samstats/star/{{sample}}_star/"
    log:
        f"{Workdir}/logs/samtools_stats/err_{{sample}}.txt"
    threads: 1
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate Samtools
        samtools stats {input.star} > {output.stats_star} 2> {log}
        plot-bamstats -p {params.out_star} {output.stats_star} 2>> {log}
        """

# ============================================================
# Filter Star Bam
# ============================================================
rule Filter_Bam:
    input:
        star=f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"
    output:
        star=temp(f"{Workdir}/alignment/star/{{sample}}.star.filter.bam")
    log:
        f"{Workdir}/logs/Filter/err_filter_{{sample}}.txt"
    threads: 1
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate Samtools
        samtools view -q 10 -b -o {output.star} {input.star} 2> {log}
        """

# ============================================================
# Sort Star Bam
# ============================================================
rule Sort_Bam:
    input:
        star=f"{Workdir}/alignment/star/{{sample}}.star.filter.bam"
    output:
        star=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
    log:
        f"{Workdir}/logs/Sort/err_sort_{{sample}}.txt"
    threads: 2
    resources:
        mem_mb=10000,
        runtime="3h"
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate Samtools
        samtools sort -o {output.star} {input.star} 2> {log}
        samtools index -b {output.star} 2>> {log}
        """






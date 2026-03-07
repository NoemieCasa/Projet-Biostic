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
        outdir = f"{Workdir}/FastQC/fastqc_raw/{{sample}}"
    log:
        f"{Workdir}/logs/fastqc_raw/err_fastqc_{{sample}}.txt"
    threads: 1
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate fastqc
        mkdir -p {params.outdir}

        # On lance FastQC
        fastqc {input.r1} {input.r2} -o {params.outdir}
        
        # Correction des noms pour correspondre au bloc output ci-dessus
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.html {output.r1_html}
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.zip {params.outdir}/{wildcards.sample}_R1_001_fastqc.zip
        
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.html {output.r2_html}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.zip {params.outdir}/{wildcards.sample}_R2_001_fastqc.zip
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
        # On utilise exactement le même chemin que l'output (sans le nom du fichier)
        outdir = f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}"
    log:
        f"{Workdir}/logs/fastqc_trimmed/err_fastqc_{{sample}}.txt"
    threads: 1
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate fastqc

        # On lance FastQC sur les deux fichiers en même temps
        fastqc {input.r1} {input.r2} -o {params.outdir} > {log} 2>&1
        """

# ============================================================
# Alignement avec Star
# ============================================================
rule Star:
    input:
        index=f"/home/iguerin2024@ec-nantes.fr/scratch/cutNtag/genome_index/star",
        r1 = f"{Workdir}/Trimming/{{sample}}_R1_001_paired.fastq.gz",
        r2 = f"{Workdir}/Trimming/{{sample}}_R2_001_paired.fastq.gz"
    output:
        bam=temp(f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"),
        summary=f"{Workdir}/alignment/star/{{sample}}.starLog.final.out"
    params:
        outpre=f"{Workdir}/alignment/star/{{sample}}.star",
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
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outpre} 2> {log}
        """

# ============================================================
# Samtools stats
# ============================================================
rule Samtools_stats:
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

# ============================================================
# fichier bam vers fichier bigwig 
# ============================================================
rule Bamtobw:
    input:
        bam=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
    output:
        bw=f"{Workdir}/bigwig/{{sample}}.bw"
    log:
        f"{Workdir}/logs/bigwig/{{sample}}.log"
    shell:
        """
		eval "$(micromamba shell hook --shell=bash)"
        micromamba activate DeepTools
        bamCoverage \
        -b {input.bam} \
        -o {output.bw} \
        --normalizeUsing CPM \
        --binSize 10 \
        > {log} 2>&1
        """

# ============================================================
# MACS2
# ============================================================
rule Macs2_callpeak:
    input:
        bam=expand(f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam", sample=SAMPLES)
    output:
        narrowpeak=f"{Workdir}/macs2/all_samples_peaks.narrowPeak"
    params:
        genome="hs"
    log:
        f"{Workdir}/logs/macs2/all_samples.log"
    shell:
        """
	    eval "$(micromamba shell hook --shell=bash)"
        micromamba activate MACS
        macs2 callpeak \
        -t {input.bam} \
        -f BAMPE \
        -g {params.genome} \
        -n all_samples \
        --outdir {Workdir}/macs2 \
        > {log} 2>&1
        """

# ============================================================
# Compute Matrix
# ============================================================
rule Compute_matrix:
    input:
        bw=expand(f"{Workdir}/bigwig/{{sample}}.bw", sample=SAMPLES),
        bed=f"{Workdir}/macs2/all_samples.bed"
    output:
        matrix=f"{Workdir}/deeptools/matrix_peaks.gz"
    shell:
        """
	micromamba activate DeepTools
        computeMatrix reference-point \
        -S {input.bw} \
        -R {input.bed} \
        -a 3000 \
        -b 3000 \
        -o {output.matrix}
        """

# ============================================================
# Plot Heatmap
# ============================================================
rule plot_heatmap:
    input:
        matrix=f"{Workdir}/deeptools/matrix_peaks.gz"
    output:
        heatmap=f"{Workdir}/deeptools/heatmap_peaks.png"
    shell:
        """
        micromamba activate DeepTools
        plotHeatmap \
        -m {input.matrix} \
        -out {output.heatmap} \
        --colorMap RdBu \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "Intensité des signaux autour des peaks"
        """
# ============================================================
# Deeptools MultiBamSummary
# ============================================================
rule multiBamSummary:
    input:
        bams=expand(f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam", sample=SAMPLES),
        bed=f"{Workdir}/macs2/all_samples.bed"
    output:
        npz=f"{Workdir}/deeptools/multiBamSummary/multibamsummary_peaks.npz",
        raw=f"{Workdir}/deeptools/multiBamSummary/multibamsummary_peaks.tab"
    log:
        f"{Workdir}/logs/deeptools/multiBamSummary_peaks.log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime="4h"
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate deeptools

        multiBamSummary BED-file \
            --BED {input.bed} \
            -b {input.bams} \
            -o {output.npz} \
            --outRawCounts {output.raw} \
            -p {threads} \
            2> {log}
        """

# ============================================================
# Deeptools plotPCA
# ============================================================
rule plotPCA:
    input:
        summary=f"{Workdir}/deeptools/multiBamSummary/multibamsummary_peaks.npz",
        bams=expand(f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam", sample=SAMPLES)
    output:
        pdf=f"{Workdir}/deeptools/PCA/PCA_peaks.pdf",
        tab=f"{Workdir}/deeptools/PCA/PCA_peaks.tab"
    log:
        f"{Workdir}/logs/deeptools/plotPCA.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime="2h"
    params:
        labels=" ".join(SAMPLES)
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate deeptools

        plotPCA \
            -in {input.summary} \
            -o {output.pdf} \
            --outFileNameData {output.tab} \
            -T "PCA of BAM samples (peaks)" \
            --labels {params.labels} \
            2> {log}
        """


##################
###MultiQC ###
##################
rule MultiQC:
	input:
		expand([f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
				f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html",
				f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R1_001_paired_fastqc.html",
				f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R2_001_paired_fastqc.html",
				f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam",
				f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam",
				f"{Workdir}/deeptools/heatmap_peaks.png"], sample=SAMPLES)
	output:
		f"{Workdir}/multiQC/multiqc_report.html"
	params:
		dir1=f"{Workdir}/FastQC/fastqc_raw/",
		dir2=f"{Workdir}/FastQC/fastqc_trimmed/",
		dir3=f"{Workdir}/alignment/star/",
		dir4=f"{Workdir}/deeptools/",
		outdir=f"{Workdir}/multiQC/"
	log:
		f"{Workdir}/logs/multiQC/err_multiQC.txt"
	shell:
		"""
		micromamba activate multiQC
		multiqc --config {Snakedir}/MultiQC_config/multiqc_config.yaml \
			--interactive \
			--outdir {params.outdir} \
			{params.dir1} {params.dir2} {params.dir3} {params.dir4} 2> {log}
		"""












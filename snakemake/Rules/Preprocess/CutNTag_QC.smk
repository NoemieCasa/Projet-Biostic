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
# Fichier bam vers fichier bigwig 
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
        narrowpeak=f"{Workdir}/macs2/all_samples_peaks.narrowPeak"
    output:
        matrix=f"{Workdir}/deeptools/matrix_peaks.gz"
    log:
        f"{Workdir}/logs/deeptools/computeMatrix.log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime="4h"
    params:
        labels=" ".join(SAMPLES)
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate DeepTools
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.narrowpeak} \
            --referencePoint center \
            -b 3000 \
            -a 3000 \
            --skipZeros \
            -p {threads} \
            -o {output.matrix} \
            2> {log}
        """

# ============================================================
# Plot Heatmap
# ============================================================
rule Plot_heatmap:
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
rule MultiBamSummary:
    input:
        bams=expand(f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam", sample=SAMPLES),
        narrowpeak=f"{Workdir}/macs2/all_samples_peaks.narrowPeak"
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
        micromamba activate DeepTools

        multiBamSummary BED-file \
            --BED {input.narrowpeak} \
            -b {input.bams} \
            -o {output.npz} \
            --outRawCounts {output.raw} \
            -p {threads} \
            2> {log}
        """
# ============================================================
# Deeptools plotPCA
# ============================================================
rule PlotPCA:
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
        micromamba activate DeepTools

        plotPCA \
            -in {input.summary} \
            -o {output.pdf} \
            --outFileNameData {output.tab} \
            -T "PCA of BAM samples (peaks)" \
            --labels {params.labels} \
            2> {log}
        """

# ============================================================
# Deeptools plotProfile
# ============================================================
rule PlotProfile:
    input:
        matrix=f"{Workdir}/deeptools/matrix_peaks.gz"
    output:
        pdf=f"{Workdir}/deeptools/plotProfile/plotProfile_peaks.pdf",
        tab=f"{Workdir}/deeptools/plotProfile/plotProfile_peaks.tab"
    log:
        f"{Workdir}/logs/deeptools/plotProfile.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime="2h"
    params:
        labels=" ".join(SAMPLES)
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate DeepTools

        plotProfile \
            -m {input.matrix} \
            -out {output.pdf} \
            --outFileNameData {output.tab} \
            --samplesLabel {params.labels} \
            -T "Average signal around peaks" \
            --perGroup \
            2> {log}
        """
# ============================================================
# MultiQC
# ============================================================
rule MultiQC:
    input:
        expand([
            f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R1_001_fastqc.html",
            f"{Workdir}/FastQC/fastqc_raw/{{sample}}/{{sample}}_R2_001_fastqc.html",
            f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R1_001_paired_fastqc.html",
            f"{Workdir}/FastQC/fastqc_trimmed/{{sample}}/{{sample}}_R2_001_paired_fastqc.html",
            f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam",
            f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam",
            f"{Workdir}/deeptools/matrix_peaks.gz",
            f"{Workdir}/deeptools/multiBamSummary/multibamsummary_peaks.npz",
            f"{Workdir}/deeptools/multiBamSummary/multibamsummary_peaks.tab",
            f"{Workdir}/deeptools/PCA/PCA_peaks.pdf",
            f"{Workdir}/deeptools/PCA/PCA_peaks.tab",
            f"{Workdir}/deeptools/plotProfile/plotProfile_peaks.pdf",
            f"{Workdir}/deeptools/plotProfile/plotProfile_peaks.tab",
            f"{Workdir}/deeptools/heatmap_peaks.png",
            f"{Workdir}/macs2/all_samples_peaks.narrowPeak",
        ], sample=SAMPLES)

    output:
        f"{Workdir}/multiQC/multiqc_report.html"

    params:
        dirs=f"{Workdir}/FastQC/fastqc_raw/ {Workdir}/FastQC/fastqc_trimmed/ {Workdir}/alignment/star/ {Workdir}/deeptools/ {Workdir}/macs2/",
        outdir=f"{Workdir}/multiQC/"

    log:
        f"{Workdir}/logs/multiQC/err_multiQC.txt"

    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate multiQC

        multiqc \
            --config {Snakedir}/MultiQC_config/multiqc_config.yaml \
            --interactive \
            --outdir {params.outdir} \
            {params.dirs} \
            2> {log}
        """

# ============================================================
# Annotation des peaks avec Homer
# ============================================================
# ============================================================
# Annotation des peaks avec Homer (Version Anti-Arobase)
# ============================================================
rule Homer_annotate_peaks:
    input:
        narrowpeak=f"{Workdir}/macs2/all_samples_peaks.narrowPeak"
    output:
        annotation=f"{Workdir}/homer/peaks_annotation.txt"
    log:
        f"{Workdir}/logs/homer/annotate_peaks.log"
    params:
        genome="hg38",
        # On enlève {sample} qui causait l'erreur
        homer_workdir="/tmp/homer_iguerin_global" 
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate homer_env

        # 1. Nettoyage et création du dossier de travail
        rm -rf {params.homer_workdir}
        mkdir -p {params.homer_workdir}
        
        # 2. Liens symboliques vers l'environnement (sans @ dans le chemin)
        ln -sfn /home/iguerin2024@ec-nantes.fr/micromamba/envs/homer_env/bin {params.homer_workdir}/bin
        ln -sfn /home/iguerin2024@ec-nantes.fr/micromamba/envs/homer_env/share/homer {params.homer_workdir}/share

        # 3. Export des variables
        export PATH="{params.homer_workdir}/bin:$PATH"
        export PERL5LIB="{params.homer_workdir}/bin:$PERL5LIB"
        export HOMER_HOME="{params.homer_workdir}/share"

        # 4. Exécution
        # On utilise le chemin relatif au bin pour tromper Perl
        {params.homer_workdir}/bin/annotatePeaks.pl \
            {input.narrowpeak} \
            {params.genome} \
            1> {output.annotation} 2> {log}

        # Nettoyage final
        rm -rf {params.homer_workdir}
        """
# ============================================================
# Extraction des catégories d'annotation pour Heatmap
# ============================================================
# ============================================================
# Extraction propre des catégories d'annotation (BED pur)
# ============================================================
rule Split_annotations_to_bed:
    input:
        annotation=f"{Workdir}/homer/peaks_annotation.txt"
    output:
        tss=f"{Workdir}/homer/split_bed/promoters.bed",
        distal=f"{Workdir}/homer/split_bed/distal_intergenic.bed"
    shell:
        """
        mkdir -p $(dirname {output.tss})

        # 1. Extraire les Promoteurs (TSS)
        # On exclut l'en-tête (PeakID), on cherche 'promoter-TSS'
        # On garde les colonnes 2, 3, 4 (Chr, Start, End) de HOMER
        awk -F'\\t' '$8 ~ /promoter-TSS/ {{print $2"\\t"$3"\\t"$4}}' {input.annotation} > {output.tss}

        # 2. Extraire le reste (Intergenic / Intron / etc.)
        # On exclut l'en-tête ET le promoter-TSS
        awk -F'\\t' '$1 !~ /PeakID/ && $8 !~ /promoter-TSS/ {{print $2"\\t"$3"\\t"$4}}' {input.annotation} > {output.distal}
        
        # Sécurité : Si un fichier est vide, deepTools va planter. 
        # On vérifie et on ajoute une ligne bidon si besoin ou on s'assure qu'ils existent.
        [ ! -s {output.tss} ] && echo "chr1\\t1\\t2" > {output.tss}
        [ ! -s {output.distal} ] && echo "chr1\\t1\\t2" > {output.distal}
        """

# ============================================================
# Matrix groupée par annotation
# ============================================================
rule Compute_matrix_annotated:
    input:
        bw=expand(f"{Workdir}/bigwig/{{sample}}.bw", sample=SAMPLES),
        regions=[f"{Workdir}/homer/split_bed/promoters.bed", f"{Workdir}/homer/split_bed/distal_intergenic.bed"]
    output:
        matrix=f"{Workdir}/deeptools/matrix_annotated.gz"
    threads: 4
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate DeepTools
        computeMatrix reference-point \
            -S {input.bw} \
            -R {input.regions} \
            --referencePoint center \
            -b 3000 -a 3000 \
            -o {output.matrix}
        """

# ============================================================
# Heatmap annotée
# ============================================================

rule Plot_heatmap_annotated:
    input:
        matrix=f"{Workdir}/deeptools/matrix_annotated.gz"
    output:
        heatmap=f"{Workdir}/deeptools/heatmap_annotated.png"
    shell:
        """
        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate DeepTools
        plotHeatmap \
            -m {input.matrix} \
            -out {output.heatmap} \
            --colorMap Viridis \
            --regionsLabel "Promoteurs" "Autres" \
            --plotTitle "Signal CutNTag par type d'annotation"
        """


















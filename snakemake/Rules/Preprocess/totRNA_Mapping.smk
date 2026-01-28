#################
###Index Fasta###
#################
rule Index_fasta:
	input:
		FastaFiles
	output:
		index=f"{Workdir}/DataBases/HG38_fasta/STAR_index/SA"
	log:
		f"{Workdir}/logs/star_index/err_star_index.txt"
	params:
		outdir=f"{Workdir}/DataBases/HG38_fasta/STAR_index/",
		gtffile=config["gtffile"]
	resources:
		mem_mb=32000,
		runtime="5h"
	threads: 
		12
	shell:
		"""
		micromamba activate STAR
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input} --sjdbGTFfile {params.gtffile} --sjdbOverhang 50 2> {log}
		"""

####################
###Star alignment###
####################	
rule star:
	input:
		index=f"{Workdir}/DataBases/HG38_fasta/STAR_index/SA",
		pair1=f"{Workdir}/trimmed_data/{{sample}}_trimmed_1P.fastq.gz",
		pair2=f"{Workdir}/trimmed_data/{{sample}}_trimmed_2P.fastq.gz"
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
	threads: 
		8
	shell:
		"""
		micromamba activate STAR
		STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.pair1} {input.pair2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outpre} 2> {log}
		"""
		
####################
###Samtools stats###
####################
rule samtools_stats:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"
	output:
		stats_star=f"{Workdir}/samstats/star/{{sample}}.star.stats"
	params:
		out_star=f"{Workdir}/samstats/star/{{sample}}_star/"
	log:
		f"{Workdir}/logs/samtools_stats/err_{{sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate Samtools
		samtools stats {input.star} > {output.stats_star} 2> {log}
		plot-bamstats -p {params.out_star} {output.stats_star} 2>> {log}
		"""

#####################
###Filter Star Bam###
#####################
rule Filter_Bam:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.starAligned.sortedByCoord.out.bam"
	output:
		star=temp(f"{Workdir}/alignment/star/{{sample}}.star.filter.bam")
	log:
		f"{Workdir}/logs/Filter/err_filter_{{sample}}.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate Samtools
		samtools view -q 10 -b -o {output.star} {input.star} 2> {log}
		"""
		
###################
###Sort Star Bam###
###################
rule Sort_Bam:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.star.filter.bam"
	output:
		star=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
	log:
		f"{Workdir}/logs/Sort/err_sort_{{sample}}.txt"
	threads: 
		2
	resources:
		mem_mb=10000,
		runtime="3h"
	shell:
		"""
		micromamba activate Samtools
		samtools sort -o {output.star} {input.star} 2> {log}
		samtools index -b {output.star} 2>> {log}
		"""




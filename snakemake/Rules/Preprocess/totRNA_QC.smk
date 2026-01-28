#################
###Make bigwig###
#################
rule Bigwig:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
	output:
		bg_star=f"{Workdir}/bigwig/star/{{sample}}.star.filter.bigwig"
	log:
		f"{Workdir}/logs/bigwig/err_bigwig_{{sample}}.txt"
	resources:
		mem_mb=10000,
		runtime="2h"
	threads: 
		1
	shell:
		"""
		micromamba activate DeepTools
		bamCoverage -b {input.star} --normalizeUsing RPKM -o {output.bg_star} 2> {log}
		"""

##################
###PlotCoverage###
##################
rule QCcoverage:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
	output:
		Raw_star=f"{Workdir}/DeepTools/star/{{sample}}/QCcoverage_{{sample}}.star.tab"
	params:
		Plot_star=f"{Workdir}/DeepTools/star/{{sample}}/QCcoverage_{{sample}}.star.png",
		Title=f"coverage_{{sample}}"
	log:
		f"{Workdir}/logs/QCcoverage/err_QCcoverage_{{sample}}.txt"
	resources:
		mem_mb=10000,
		runtime="2h"
	threads: 
		1
	shell:
		"""
		micromamba activate DeepTools
		plotCoverage -b {input.star} --skipZeros --centerReads --plotFile {params.Plot_star} --plotTitle {params.Title} --outRawCounts {output.Raw_star} 2> {log}
		"""

#######################
###BamPEfragmentSize###
#######################
rule BamPEFRagmentsize:
	input:
		star=f"{Workdir}/alignment/star/{{sample}}.star.filter.sort.bam"
	output:
		table_star=f"{Workdir}/DeepTools/star/{{sample}}/fragmentSize_Table.{{sample}}.star.txt",
		rawLen_star=f"{Workdir}/DeepTools/star/{{sample}}/fragmentSize_RawLen.{{sample}}.star.txt"
	params:
		Hist_star=f"{Workdir}/DeepTools/star/{{sample}}/fragmentSize_{{sample}}.star.png",
		Title=f"fragmentSize_{{sample}}",
		sample=f"{{sample}}"
	log:
		f"{Workdir}/logs/PEFragments/err_PEFragments_{{sample}}.txt"
	resources:
		mem_mb=10000,
		runtime="2h"
	threads: 
		1
	shell:
		"""
		micromamba activate DeepTools
		bamPEFragmentSize -hist {params.Hist_star} --maxFragmentLength 1000 -T {params.Title} --table {output.table_star} --outRawFragmentLengths {output.rawLen_star} -b {input.star} --samplesLabel {params.sample} 2> {log}
		"""
		
##################
###MultiQC star###
##################
rule multiqc_star:
	input:
		expand([f"{Workdir}/FastQC/trimmed_data/{{sample}}_trimmed_1P_fastqc.zip",
				f"{Workdir}/FastQC/trimmed_data/{{sample}}_trimmed_2P_fastqc.zip",
				f"{Workdir}/alignment/star/{{sample}}.starLog.final.out",
				f"{Workdir}/samstats/star/{{sample}}.star.stats",
				f"{Workdir}/DeepTools/star/{{sample}}/QCcoverage_{{sample}}.star.tab",
				f"{Workdir}/DeepTools/star/{{sample}}/fragmentSize_RawLen.{{sample}}.star.txt",
				f"{Workdir}/DeepTools/star/{{sample}}/fragmentSize_Table.{{sample}}.star.txt",
				f"{Workdir}/featureCounts/star/{{sample}}.star_rawcounts.txt.summary"],				
				sample=SAMPLES)
	output:
		f"{Workdir}/multiQC/star/multiqc_report.html"
	params:
		dir1=f"{Workdir}/FastQC/trimmed_data/",
		dir2=f"{Workdir}/alignment/star/",
		dir3=f"{Workdir}/samstats/star/",
		dir4=f"{Workdir}/DeepTools/star/",
		dir5=f"{Workdir}/featureCounts/star/",
		outdir=f"{Workdir}/multiQC/star/"
	log:
		f"{Workdir}/logs/multiQC_star/err_multiQC.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate multiQC
		multiqc --config {Snakedir}/MultiQC_config/multiqc_config.yaml \
			--interactive \
			--outdir {params.outdir} \
			{params.dir1} {params.dir2} {params.dir3} {params.dir4} {params.dir5} 2> {log}
		"""

#########################
###MultiQC on raw data###
#########################
rule multiQC_raw:
	input:
		FilesForMultiQC
	output:
		f"{Workdir}/multiQC/raw/multiqc_report.html"
	params:
		dir1=f"{Workdir}/FastQC/raw_data/",
		outdir=f"{Workdir}/multiQC/raw/"
	log:
		f"{Workdir}/logs/multiQC_raw/err_multiQC.txt"
	threads: 
		1
	shell:
		"""
		micromamba activate multiQC
		multiqc --config {Snakedir}/MultiQC_config/multiqc_config.yaml \
			--interactive \
			--outdir {params.outdir} \
			{params.dir1} 2> {log}
		"""

#################
###Make bigwig###
#################
rule Bigwig_all:
	input:
		Treat=TreatmentFiles_B,
		Cont=ControlFiles_B,
		Ref=RefFiles_B
	output:
		bw_Treat=f"{Workdir}/bigwig/Rejection.star.bigwig",
		bw_Cont=f"{Workdir}/bigwig/NoRejection.star.bigwig",
		bw_Ref=f"{Workdir}/bigwig/Reference.star.bigwig"
	params:
		bg_Treat=f"{Workdir}/bigwig/Rejection.star.bedgraph",
		bg_Cont=f"{Workdir}/bigwig/NoRejection.star.bedgraph",
		bg_Ref=f"{Workdir}/bigwig/Reference.star.bedgraph",
		bg_Treat_sorted=f"{Workdir}/bigwig/Rejection.star.sorted.bedgraph",
		bg_Cont_sorted=f"{Workdir}/bigwig/NoRejection.star.sorted.bedgraph",
		bg_Ref_sorted=f"{Workdir}/bigwig/Reference.star.sorted.bedgraph",
		Chrom_Sizes=config["chrom_sizes"]
	log:
		f"{Workdir}/logs/bigwig/err_bigwig_all.txt"
	resources:
		mem_mb=16000,
		runtime="2h"
	threads: 
		1
	shell:
		"""
		micromamba activate Bigwigmerge
		bigWigMerge {input.Treat} {params.bg_Treat} 2> {log}
		sort -k1,1 -k2,2n {params.bg_Treat} > {params.bg_Treat_sorted} 2>> {log}
		bedGraphToBigWig {params.bg_Treat_sorted} {params.Chrom_Sizes} {output.bw_Treat} 2>> {log}
		bigWigMerge {input.Cont} {params.bg_Cont} 2>> {log}
		sort -k1,1 -k2,2n {params.bg_Cont} > {params.bg_Cont_sorted} 2>> {log}
		bedGraphToBigWig {params.bg_Cont_sorted} {params.Chrom_Sizes} {output.bw_Cont} 2>> {log}
		bigWigMerge {input.Ref} {params.bg_Ref} 2>> {log}
		sort -k1,1 -k2,2n {params.bg_Ref} > {params.bg_Ref_sorted} 2>> {log}
		bedGraphToBigWig {params.bg_Ref_sorted} {params.Chrom_Sizes} {output.bw_Ref} 2>> {log}
		"""
		
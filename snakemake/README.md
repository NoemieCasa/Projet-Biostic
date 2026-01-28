# IMoKiT_totRNA-seq

## Getting started


- [ ] Download refseq databases for gene annotations at [Refseq](https://www.ncbi.nlm.nih.gov/refseq/).
- [ ] Download the human genome data  in fasta format and chromosome size data on [ucsc](https://hgdownload.soe.ucsc.edu/downloads.html).
- [ ] Download the FANTOM5 database for enhancer expressions [FANTOM5](https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/)


- [ ] You will need micromamba to use snakemake (see how to install it on glicid [micromamba](https://doc.glicid.fr/GLiCID-PUBLIC/environment/micromamba.html)).
- [ ] To manage module conflicts, you will need to install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). (only bird users)


- [ ] Clone [trinity](https://github.com/trinityrnaseq/trinityrnaseq) github and set trinity path in config.json

### Singularity images
- [ ] Download singularity image for trinity (change the path in config.json) [trinity](https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/)
- [ ] Metascape for gsea
- [ ] CibersortX for deonvolution

- [ ] Create all new conda environments form /Conda to install snakemake (not definitive, micromamba is not included in snakemake for now) -> thinking about switching to singularity
```
micromamba env create -f Conda/IMoKiT.yml
```

- [ ] Change path and option for your working or raw directory, sample names etc in Config.json
- [ ] Create logs and logs_gen directories (do not forget to change the path)
- [ ] PATH to metascape (+ how to install it)

## Run Snakemake

Run snakemake
```
conda activate IMoKiT
snakemake --profile profile/ > logs_gen/P01-01_Snakemake.o.txt 2> logs_gen/P01-01_Snakemake.e.txt
```

## Visualization

We removed multiQC and all rules to have a better vizualisation

```
snakemake --dag --use-conda | 
	grep -v "\-> 0\|0\[label = \"all\"" | 
	grep -v "\-> 16\|16\[label = \"multiqc_star\"" |
	grep -v "\-> 25\|25\[label = \"multiQC_raw\"" |
	grep -v "26 \->\|26\[label = \"Fastqc_raw" |
	grep -v "27 \->\|27\[label = \"Fastqc_raw" |
	dot -Tsvg > dag.svg
```

## Log
- 2025-09-04 Initial commit for preprocess raw data
- 2025-09-08 Change the q for priority and account to team3
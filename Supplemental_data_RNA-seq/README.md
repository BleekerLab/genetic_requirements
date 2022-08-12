# RNA-seq analysis

This folder contains the supplemental RNA-seq processed data used for the differential expression analysis and other representations (PCA, heatmaps) for this paper. 

<!-- MarkdownTOC autolink="true" levels="1,2" -->

- [Reports from the analysis](#reports-from-the-analysis)
	- [MultiQC report from the individual fastp reports](#multiqc-report-from-the-individual-fastp-reports)
	- [Mapping report](#mapping-report)
- [Code and data used](#code-and-data-used)
	- [Sequencing files](#sequencing-files)
	- [Software](#software)
	- [Genomic references](#genomic-references)

<!-- /MarkdownTOC -->


# Reports from the analysis

## MultiQC report from the individual fastp reports

Open the [MultiQC report](./multiqc_report.html) in a web browser. 

## Mapping report

See [mapping summary table](./mapping_summary.csv)

# Code and data used

## Sequencing files
- F1, Elite line and S. habrochaites PI127826 mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.3603229) 
- F2 lines mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.4491747)

The F2-28 replicates 1 to 4 were pooled together to constitute one unique F2-28 individual plant.  
The md5 checksum of the pooled F2-28 sample is `74032ef1726a9013bf4b53f7475fbc2f`. 

## Software 
Snakemake mRNA-seq pipeline v0.3.8 based on STAR: [link](https://github.com/BleekerLab/snakemake_rnaseq/releases/tag/v0.3.8)


## Genomic references
The tomato annotation ITAG2.4 was used. Specifally the cDNA fasta file was used.   

- [SOL Genomics link](https://solgenomics.net/organism/Solanum_lycopersicum/genome). Go to "FTP Site" then "annotation" and "ITAG2.4" for the "ITAG2.4_gene_models.gff3".     
- [Custom Zenodo record for the tomato genomic references](https://zenodo.org/record/5578601). The genome is "S_lycopersicum_chromosomes.2.50.fa.gz" and the annotation "ITAG2.4_gene_models.gff".  

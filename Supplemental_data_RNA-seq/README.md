# RNA-seq analysis

This folder contains the `raw_counts.tsv` and the `scaled_counts.tsv` that are used for the differential expression analysis and other representations (PCA, heatmaps) for this paper. 

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

See [mapping summary table](./mapping_results.csv)

| n_targets 	| n_bootstraps 	| n_processed 	| n_pseudoaligned 	| n_unique 	| p_pseudoaligned 	| p_unique 	| sample      	|
|-----------	|--------------	|-------------	|-----------------	|----------	|-----------------	|----------	|-------------	|
| 34075     	| 100          	| 47839301    	| 43861023        	| 43085770 	| 91.7            	| 90.1     	| Elite_02    	|
| 34075     	| 100          	| 25763854    	| 22315511        	| 22086331 	| 86.6            	| 85.7     	| F2-127      	|
| 34075     	| 100          	| 23387126    	| 20474496        	| 19872926 	| 87.5            	| 85       	| F2-151      	|
| 34075     	| 100          	| 83670152    	| 57901546        	| 55706238 	| 69.2            	| 66.6     	| F2-28       	|
| 34075     	| 100          	| 29535331    	| 25972873        	| 25559982 	| 87.9            	| 86.5     	| F2-411      	|
| 34075     	| 100          	| 27980463    	| 24514301        	| 24165926 	| 87.6            	| 86.4     	| F2-445      	|
| 34075     	| 100          	| 28190078    	| 24914196        	| 24645285 	| 88.4            	| 87.4     	| F2-73       	|
| 34075     	| 100          	| 48627207    	| 41552140        	| 40707092 	| 85.5            	| 83.7     	| PI127826    	|
| 34075     	| 100          	| 32277727    	| 28364133        	| 27987404 	| 87.9            	| 86.7     	| PI127826_F1 	|

# Code and data used

## Sequencing files
- F1, Elite line and S. habrochaites PI127826 mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.3603229) 
- F2 lines mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.4491747)

The F2-28 replicates 1 to 4 were pooled together to constitute one unique F2-28 individual plant.  
The md5 checksum of the pooled F2-28 sample is `74032ef1726a9013bf4b53f7475fbc2f`. 

## Software 
Snakemake mRNA-seq pipeline based on Kallisto: [link](https://zenodo.org/record/4091104)


## Genomic references
The tomato annotation ITAG4.0 was used. Specifally the cDNA fasta file was used.   

[SOL Genomics link](ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_cDNA.fasta)  
[Custom Zenodo record](https://10.5281/zenodo.4321000)

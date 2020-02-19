# Figure 6

## Data provenance
The `abundance_tidy.tsv` file contains the transcript scaled counts based on mRNA-seq files (see below) and on the ITAG4.0 transcriptome of S. lycopersicum Heinz1706.

mRNA-seq fastq files archives:
* 2017-11-23: _S. lycopersicum_ and _S. habrochaites_ parental lines and F1s: https://zenodo.org/record/3610267  
* 2019-03-25: selected F2 "active" and "lazy" _S. lycopersicum_ Elite x _S. habrochaites_ PI127826 F2 lines: https://zenodo.org/record/3610279

These mRNA-seq fastq files were processed using a Snakemake pipeline available here: https://github.com/BleekerLab/rnaseq-analysis-kallisto-sleuth/releases/tag/v0.2.2  

The Solanum lycopersicum ITAG4.0 transcriptome `ITAG4.0_cDNA.fasta` was used as the reference and can be retrieved here: https://zenodo.org/record/3613695   
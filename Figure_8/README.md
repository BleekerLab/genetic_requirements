# Differential expression

![Figure 8](./figure8.png)


# Data provenance

## Derived data: count tables
The `raw_counts.txt` are derived from mRNA-seq fastq files (see "mRNA-seq fastq files" below) using a dedicated pipeline
(see "Software" section below). 
These files contain raw counts based on _S. lycopersicum_ Heinz1706 genome ITAG4.0.

The `raw_counts.txt` is available on [Zenodo](https://doi.org/10.5281/zenodo.3959143).

## Raw data: mRNA-seq fastq files

mRNA-seq fastq files archives:
* _S. lycopersicum_ and _S. habrochaites_ parental lines and F1s: [Link to Zenodo archive](https://zenodo.org/record/3610267).  
* Selected F2 "active" and "lazy" lines from a cross between _S. lycopersicum_ Elite x _S. habrochaites_ PI127826: [Link to Zenodo archive](https://zenodo.org/record/3610279).

## Software
A pipeline to go from RNA-seq fastq files to counts was used and can be found on [GitHub](https://github.com/BleekerLab/snakemake_rnaseq_to_counts/releases/tag/v0.3.0).


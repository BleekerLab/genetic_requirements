# Differential expression analysis of active and lazy genotypes

The description of the fastq mRNA-seq to Kallisto pseudoalignments is described in the [`Supplemental_data_RNA-seq/` folder.](../Supplemental_data_RNA-seq/)  

# Chromomaps

The [chromoMap](https://cran.r-project.org/web/packages/chromoMap/index.html) package depends on R >4.0. Since I still use R 3.6 for compatibility reasons, I will use a Docker image to create the chromomaps. 

## RStudio 4.1 docker image
To create the chromomaps figure, first get the RStudio docker image for version 4.1: `docker pull rocker/rstudio:4.1`.

Inside `Figure_08_diff_expr/`, do:  
```bash
docker run --detach -p 8787:8787 -e DISABLE_AUTH=true -v $PWD:/home/rstudio/workspace/ rocker/rstudio:4.1
```
Then go to [http://localhost:8787/](http://localhost:8787/). 

## Install chromoMap

1. `install.packages("chromoMap")`
2. 

## Chromosome sizes

These were obtained using samtools:  
```bash
samtools faidx S_lycopersicum_chromosomes.2.50.fa
cut -f1,2 S_lycopersicum_chromosomes.2.50.fa.fai > sizes.genome
```
Followed by:
```bash
awk '{print $1 "\t" 1 "\t" $2}' sizes.genome 
```
This gives:  
```
SL2.50ch00	1	21805821
SL2.50ch01	1	98543444
SL2.50ch02	1	55340444
SL2.50ch03	1	70787664
SL2.50ch04	1	66470942
SL2.50ch05	1	65875088
SL2.50ch06	1	49751636
SL2.50ch07	1	68045021
SL2.50ch08	1	65866657
SL2.50ch09	1	72482091
SL2.50ch10	1	65527505
SL2.50ch11	1	56302525
SL2.50ch12	1	67145203
```

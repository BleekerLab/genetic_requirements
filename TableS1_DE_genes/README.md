# Table S1
Differentially expressed genes between the "lazy" and "active" trichome types. Here is the experimental design file:  

| sample | trichome | path                              |
|--------|----------|-----------------------------------|
| F2-127 | active   | TableS1_DE_genes/kallisto/F2-127/ |
| F2-73  | active   | TableS1_DE_genes/kallisto/F2-73/  |
| F2-28  | active   | TableS1_DE_genes/kallisto/F2-28/  |
| F2-411 | .lazy    | TableS1_DE_genes/kallisto/F2-411/ |
| F2-445 | .lazy    | TableS1_DE_genes/kallisto/F2-445/ |
| F2-151 | .lazy    | TableS1_DE_genes/kallisto/F2-151/ |

## Graphical visualisation

### PCA
![PCA score plot]("./plot_pca.png")
![PCA loading plot]("./plot_loadings.png")

### Volcano
![Volcano plot]("./plot_volcano.png")



## Method
* Sleuth differential expression method: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
* Using the F2 lines phenotype "lazy" or "active" as the explanatory variable.  
* Annotations are from ITAG4.1 (Sol genomics)

## Results
The list of differentially expressed genes with their annotation is available in `differentials.tsv`. 
Column explanation are:  
* target_id: transcript name, e.g. "ENST#####" (dependent on the transcriptome used in kallisto). If gene_mode is TRUE, this will instead be the IDs specified by the obj$gene_column from obj$target_mapping.
* pval: p-value of the chosen model
* qval: false discovery rate adjusted p-value, using Benjamini-Hochberg (see p.adjust)
* b: 'beta' value (effect size). Technically a biased estimator of the fold change. Only seen with Wald test results.
* se_b: standard error of the beta. Only seen with Wald test results.
* mean_obs: mean of natural log counts of observations
* var_obs: variance of observation
* tech_var: technical variance of observation from the bootstraps (named 'sigma_q_sq' if rename_cols is FALSE)
* sigma_sq: raw estimator of the variance once the technical variance has been removed
* smooth_sigma_sq: smooth regression fit for the shrinkage estimation
* final_simga_sq: max(sigma_sq, smooth_sigma_sq); used for covariance estimation of beta (named 'smooth_sigma_sq_pmax' if rename_cols is FALSE)


## Data provenance
The `abundance_tidy.tsv` file contains the transcript scaled counts based on mRNA-seq files (see below) and on the ITAG4.0 transcriptome of S. lycopersicum Heinz1706.

mRNA-seq fastq files archives:
* 2017-11-23: _S. lycopersicum_ and _S. habrochaites_ parental lines and F1s: https://zenodo.org/record/3610267  
* 2019-03-25: selected F2 "active" and "lazy" _S. lycopersicum_ Elite x _S. habrochaites_ PI127826 F2 lines: https://zenodo.org/record/3610279

These mRNA-seq fastq files were processed using a Snakemake pipeline available on [GitHub](https://github.com/BleekerLab/rnaseq-analysis-kallisto-sleuth/releases/tag/v0.2.2) and [Zenodo](https://doi.org/10.5281/zenodo.3627098).  

The Solanum lycopersicum ITAG4.0 transcriptome `ITAG4.0_cDNA.fasta` was used as the reference and can be retrieved here: https://zenodo.org/record/3613695   



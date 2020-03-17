# iPath v3.0 analysis 
The KEGG KO identifiers (KEGG Orthology id) were retrieved for the 
134 differential genes (qval < 0.01) found after the
[DESeq2 analysis](https://github.com/BleekerLab/genetic_requirements/tree/master/TableS1_DE_genes). 

The [BlastKOALA tool](https://www.kegg.jp/blastkoala/) was used to retrieve their KO identifiers. Then, these identifiers were submitted to the [iPath v3.0 tool](https://pathways.embl.de/) to create metabolic maps. 

## Shared metabolic maps
All created maps can be retrieved here: https://pathways.embl.de/shared/mgalland

## Data used
1. For BlastKOALA, the `differential_genes_protein_sequences.fasta` file was used (ITAG4.0 proteins).
2. For iPath v3.0, the `kegg` column from the `differentially_expressed_active_vs_lazy_qval001.tsv` was used as it contains the KO identifiers.
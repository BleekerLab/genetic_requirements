# RNA-seq analysis

This folder contains the `raw_counts.tsv` and the `scaled_counts.tsv` that are used for the differential expression analysis and other representations (PCA, heatmaps) for this paper. 


# Reports from the analysis

## FastQC metrics

## Mapping report

|                                                      | Elite_02   | PI127826_F1 | PI127826   | F2-127     | F2-151     | F2-28      | F2-411     | F2-445     | F2-73      |
|------------------------------------------------------|------------|-------------|------------|------------|------------|------------|------------|------------|------------|
|                                  Started job   on| | ########## | ##########  | ########## | ########## | ########## | ########## | ########## | ########## | ########## |
|                              Started mapping   on| | ########## | ##########  | ########## | ########## | ########## | ########## | ########## | ########## | ########## |
|                                     Finished   on| | ########## | ##########  | ########## | ########## | ########## | ########## | ########## | ########## | ########## |
|          Mapping speed, Million of reads per hour| | 167.58     | 202.43      | 136.71     | 86.31      | 82.39      | 79.96      | 79.51      | 90.05      | 83.31      |
|                           Number of input   reads| | 46644373   | 31377108    | 46672187   | 24479621   | 22153589   | 29009113   | 28113938   | 26666261   | 26845756   |
|                       Average input read   length| | 199        | 199         | 199        | 301        | 301        | 301        | 301        | 301        | 301        |
|                                     UNIQUE   READS:  |            |             |            |            |            |            |            |            |            |
|                    Uniquely mapped reads   number| | 42782744   | 24236892    | 28265618   | 13857720   | 14070513   | 18958552   | 16823157   | 15967434   | 16043053   |
|                         Uniquely mapped reads   %| | 91.72%     | 77.24%      | 60.56%     | 56.61%     | 63.51%     | 65.35%     | 59.84%     | 59.88%     | 59.76%     |
|                           Average mapped   length| | 199.78     | 199.72      | 199.64     | 301.23     | 301.34     | 301.34     | 301.34     | 301.40     | 301.30     |
|                        Number of splices:   Total| | 28104611   | 17655951    | 21320568   | 15319703   | 14830673   | 20535886   | 17892065   | 16715157   | 17542735   |
|             Number of splices: Annotated   (sjdb)| | 26289885   | 16580440    | 20196897   | 14646966   | 14105508   | 19657881   | 17074535   | 16004176   | 16796868   |
|                        Number of splices:   GT/AG| | 27775241   | 17431312    | 21047224   | 15126569   | 14608544   | 20269859   | 17642438   | 16495105   | 17313149   |
|                        Number of splices:   GC/AG| | 252232     | 169590      | 188272     | 153394     | 182744     | 218191     | 203535     | 175063     | 186223     |
|                        Number of splices:   AT/AC| | 9931       | 8243        | 8075       | 5990       | 6917       | 7342       | 6632       | 6517       | 6454       |
|                Number of splices:   Non-canonical| | 67207      | 46806       | 76997      | 33750      | 32468      | 40494      | 39460      | 38472      | 36909      |
|                       Mismatch rate per base,   %| | 0.17%      | 0.52%       | 1.00%      | 0.56%      | 0.43%      | 0.44%      | 0.50%      | 0.50%      | 0.51%      |
|                          Deletion rate per   base| | 0.01%      | 0.04%       | 0.08%      | 0.05%      | 0.03%      | 0.03%      | 0.04%      | 0.04%      | 0.04%      |
|                         Deletion average   length| | 1.55       | 2.35        | 2.73       | 2.59       | 2.36       | 2.53       | 2.50       | 2.42       | 2.69       |
|                         Insertion rate per   base| | 0.02%      | 0.04%       | 0.07%      | 0.04%      | 0.03%      | 0.03%      | 0.04%      | 0.04%      | 0.04%      |
|                        Insertion average   length| | 1.32       | 2.22        | 3.22       | 2.98       | 2.42       | 2.75       | 2.57       | 2.40       | 3.19       |
|                              MULTI-MAPPING   READS:  |            |             |            |            |            |            |            |            |            |
|           Number of reads mapped to multiple loci| | 1477766    | 672097      | 897904     | 388177     | 766710     | 636434     | 772267     | 1185883    | 636531     |
|              % of reads mapped to multiple   loci| | 3.17%      | 2.14%       | 1.92%      | 1.59%      | 3.46%      | 2.19%      | 2.75%      | 4.45%      | 2.37%      |
|           Number of reads mapped to too many loci| | 7933       | 3899        | 3929       | 2501       | 2360       | 10057      | 7225       | 2767       | 10192      |
|              % of reads mapped to too many   loci| | 0.02%      | 0.01%       | 0.01%      | 0.01%      | 0.01%      | 0.03%      | 0.03%      | 0.01%      | 0.04%      |
|                                   UNMAPPED   READS:  |            |             |            |            |            |            |            |            |            |
|          % of reads unmapped: too many mismatches| | 2.20%      | 15.74%      | 31.40%     | 35.66%     | 27.71%     | 28.17%     | 32.09%     | 29.52%     | 32.98%     |
|                  % of reads unmapped: too   short| | 0.41%      | 2.71%       | 4.60%      | 4.42%      | 3.89%      | 3.19%      | 3.61%      | 3.33%      | 3.61%      |
|                      % of reads unmapped:   other| | 2.48%      | 2.15%       | 1.50%      | 1.72%      | 1.41%      | 1.05%      | 1.69%      | 2.82%      | 1.23%      |
|                                   CHIMERIC   READS:  |            |             |            |            |            |            |            |            |            |
|                        Number of chimeric   reads| | 0          | 0           | 0          | 0          | 0          | 0          | 0          | 0          | 0          |
|                             % of chimeric   reads| | 0.00%      | 0.00%       | 0.00%      | 0.00%      | 0.00%      | 0.00%      | 0.00%      | 0.00%      | 0.00%      |

# Code and data used

## Sequencing files
- F1, Elite line and S. habrochaites PI127826 mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.3603229) 
- F2 lines mRNA-seq fastq files from stem trichomes: [link](https://doi.org/10.5281/zenodo.3610278)

## Software 
Snakemake mRNA-seq pipeline: [link](https://zenodo.org/record/4034215)

## Genomic references
Tomato assembly and annotation SL4.0


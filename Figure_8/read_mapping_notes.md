# Mapping rates notes

No mismatch allowed: a read has to match perfectly to the genome reference.  
No multimapper allowed: the read has to have one unique matching position in the genome. 

# Solanum lycopersicum Heinz1706 genome reference

## Elite line parent 
Number of input reads:	46,644,373    
Average input read length:	199    
*UNIQUE READS*    
- Uniquely mapped reads number:	32,700,930
- Uniquely mapped reads %:	**70.11%**  
*MULTI-MAPPING READS*    
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	1121965
- % of reads mapped to too many loci:	2.41%  
*UNMAPPED READS*  
- % of reads unmapped: too many mismatches:	24.37%
- % of reads unmapped: too short:	0.64%
- % of reads unmapped: other:	2.48%

## PI127826 parent 
Number of input reads:	46,672,187  
Average input read length:	199  
*UNIQUE READS*   
- Uniquely mapped reads number:	4,413,619
- Uniquely mapped reads %:	**9.46%**  
*MULTI-MAPPING READS*    
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	115617
- % of reads mapped to too many loci:	0.25%   
*UNMAPPED READS*  
- % of reads unmapped: too many mismatches:	**81.24%**  
- % of reads unmapped: too short:	7.55%
- % of reads unmapped: other:	1.50%

## F1 Elite x PI127826
Number of input reads:	31,377,108
Average input read length:	199  
*UNIQUE READS*   
- Uniquely mapped reads number:	12,377,797
- Uniquely mapped reads %:	**39.45%**    
*MULTI-MAPPING READS*   
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	332767
- % of reads mapped to too many loci:	1.06%  
*UNMAPPED READS*    
- % of reads unmapped: too many mismatches:	*53.04%*
- % of reads unmapped: too short:	4.30%
- % of reads unmapped: other:	2.15%


# Solanum habrochaites PI127826 genome reference

## Elite line parent
Number of input reads:	46,644,373
Average input read length:	199
*UNIQUE READS*    
- Uniquely mapped reads number:	3643260
- Uniquely mapped reads %:	**7.81%**    
*MULTI-MAPPING READS*   
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	1055218
- % of reads mapped to too many loci:	2.26%  
*UNMAPPED READS*  
- % of reads unmapped: too many mismatches:	**74.13%**  
- % of reads unmapped: too short:	15.61%
- % of reads unmapped: other:	0.19%

## PI127826 parent
Number of input reads:	46,672,187
Average input read length:	199  
*UNIQUE READS*
- Uniquely mapped reads number:	26212590
- Uniquely mapped reads %:	**56.16%**  
*MULTI-MAPPING READS*
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	1433662
- % of reads mapped to too many loci:	3.07%
*UNMAPPED READS*
- % of reads unmapped: too many mismatches:	**35.63%** (segregating PI127826 pop?) 
- % of reads unmapped: too short:	5.06%
- % of reads unmapped: other:	0.08%

## F1 Elite x PI127826
Number of input reads:	31,377,108
Average input read length:	199
*UNIQUE READS*  
- Uniquely mapped reads number:	9,848,760
- Uniquely mapped reads %:	**31.39%**    
*MULTI-MAPPING READS*  
- Number of reads mapped to multiple loci:	0
- % of reads mapped to multiple loci:	0.00%
- Number of reads mapped to too many loci:	625003
- % of reads mapped to too many loci:	1.99%  
*UNMAPPED READS*  
- % of reads unmapped: too many mismatches:	**56.25%**    
- % of reads unmapped: too short:	10.22%
- % of reads unmapped: other:	0.15%

# Other things to keep
Trick to remove leading whitespaces from log files.
```
sed "s/^[ \t]*//" Elite_02_Log.final.out > Elite_02_Log.final.parsed.out 
```


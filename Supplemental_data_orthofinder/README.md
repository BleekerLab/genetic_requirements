# Finding orthologues of MEP and MVA pathway genes


## Orthofinder

Orthofinder v2.5.4 was used to retrieve the orthologues from S. lycopersicum Heinz1706 and S. habrochaites PI127826. 
Default parameters were used. 

From all results, the `ITAG4.0_proteins__v__Solanum_habrochaites_PI127826_protein_Dovetails_2021.tsv` result file was used that contains the orthologues of _S. lycopersicum_ proteins in _S. habrochaites_ PI127826.  

### Datasets

* `ITAG4.0_proteins.fasta` for lycopersicum   
* `Solanum_habrochaites_PI127826_protein_Dovetails_2021.fasta` for habrochaites

Zenodo doi: [https://zenodo.org/record/5528738](https://zenodo.org/record/5528738) 

### Bash script to retrieve the orthologues given a list of ids

```bash
bash retrieve_habro_orthologues.sh
```

## MEP/MVA Orthologues

| Orthogroup 	| ITAG4.0_protein_id 	| PI127826_protein_id 	| locus          	| function 	| pathway 	|
|------------	|--------------------	|---------------------	|----------------	|----------	|---------	|
| OG0006942  	| Solyc11g010850.2.1 	| ANN32664-RA         	| Solyc11g010850 	| DXS      	| MEP     	|
| OG0015177  	| Solyc08g066950.3.1 	| ANN37024-RA         	| Solyc08g066950 	| DXS      	| MEP     	|
| OG0010934  	| Solyc03g114340.3.1 	| ANN07077-RA         	| Solyc03g114340 	| DXR      	| MEP     	|
| OG0008607  	| Solyc01g102820.4.1 	| ANN00805-RA         	| Solyc01g102820 	| MCT      	| MEP     	|
| OG0007696  	| Solyc01g009010.3.1 	| ANN04231-RA         	| Solyc01g009010 	| CMK      	| MEP     	|
| OG0015573  	| Solyc08g081570.3.1 	| ANN28001-RA         	| Solyc08g081570 	| MDS      	| MEP     	|
| OG0017739  	| Solyc11g069380.2.1 	| ANN30757-RA         	| Solyc11g069380 	| HDS      	| MEP     	|
| OG0008890  	| Solyc01g109300.3.1 	| ANN00002-RA         	| Solyc01g109300 	| HDR      	| MEP     	|
| OG0002050  	| Solyc08g005677.1.1 	| ANN30391-RA         	| Solyc08g005677 	| TPS20    	| MEP     	|
| OG0002050  	| Solyc08g005670.2.1 	| ANN30391-RA         	| Solyc08g005670 	| TPS20    	| MEP     	|
| OG0002050  	| Solyc08g005640.4.1 	| ANN30391-RA         	| Solyc08g005640 	| TPS20    	| MEP     	|
| OG0002713  	| Solyc05g017760.4.1 	| ANN20601-RA         	| Solyc05g017760 	| AACT     	| MVA     	|
| OG0002713  	| Solyc07g045350.4.1 	| ANN14752-RA         	| Solyc07g045350 	| AACT     	| MVA     	|
| OG0001070  	| Solyc08g007790.3.1 	| ANN30258-RA         	| Solyc08g007790 	| HMGS     	| MVA     	|
| OG0001070  	| Solyc08g080170.3.1 	| ANN28112-RA         	| Solyc08g080170 	| HMGS     	| MVA     	|
| OG0001070  	| Solyc12g056450.2.1 	| ANN33562-RA         	| Solyc12g056450 	| HMGS     	| MVA     	|
| OG0001184  	| Solyc03g032020.4.1 	| ANN12039-RA         	| Solyc03g032020 	| HMGR     	| MVA     	|
| OG0001184  	| Solyc02g038740.4.1 	| ANN12039-RA         	| Solyc02g038740 	| HMGR     	| MVA     	|
| OG0001184  	| Solyc02g082260.3.1 	| ANN12039-RA         	| Solyc02g082260 	| HMGR     	| MVA     	|
| OG0008473  	| Solyc01g098840.3.1 	| ANN01059-RA         	| Solyc01g098840 	| MVK      	| MVA     	|
| OG0002806  	| Solyc08g076140.3.1 	| ANN28420-RA         	| Solyc08g076140 	| pMVK     	| MVA     	|
| OG0002806  	| Solyc06g066310.3.1 	| ANN18532-RA         	| Solyc06g066310 	| pMVK     	| MVA     	|
| OG0003756  	| Solyc11g007020.2.1 	| ANN32860-RA         	| Solyc11g007020 	| MDC      	| MVA     	|
| OG0003756  	| Solyc04g009650.4.1 	| ANN32860-RA         	| Solyc04g009650 	| MDC      	| MVA     	|
| OG0012979  	| Solyc05g055760.3.1 	| ANN22251-RA         	| Solyc05g055760 	| IDI2     	| MVA     	|
| OG0015303  	| Solyc08g075390.4.1 	| ANN28488-RA         	| Solyc08g075390 	| Nudix    	| MVA     	|
| OG0003739  	| Solyc04g005530.3.1 	| ANN25256-RA         	| Solyc04g005530 	| IPK      	| MVA     	|
| OG0003739  	| Solyc04g005520.3.1 	| ANN25256-RA         	| Solyc04g005520 	| IPK      	| MVA     	|
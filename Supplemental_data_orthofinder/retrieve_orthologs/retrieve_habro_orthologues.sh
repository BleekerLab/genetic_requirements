#!/bin/bash

FILENAME="ids_precursor_genes.tsv"

IDS=$(cat $FILENAME)

for ID in $IDS;
do 
	grep "$ID" ITAG4.0_proteins__v__Solanum_habrochaites_PI127826_protein_Dovetails_2021.tsv
done 


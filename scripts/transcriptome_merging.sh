#! /bin/bash

## Authors: Francisco J. Romero-Campero
##          Ana Belen Romero-Losada
## Contact: Francisco J. Romero-campero - fran@us.es

EXP_FOLDER=$1
ANNOTATION=$2

## Access results folder
cd ${EXP_FOLDER}/results

## Merging sample transcriptomes
stringtie --merge -G $ANNOTATION -o stringtie_merged.gtf merge_list.txt

## Comparing our assembly with the reference
gffcompare -r $ANNOTATION -G -o comparison stringtie_merged.gtf

## Generation of read count and TPM
prepDE.py -i ../samples

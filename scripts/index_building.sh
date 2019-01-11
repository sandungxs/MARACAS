## Author: Francisco J. Romero-Campero
## Date: January 2019
## Email: fran@us.es

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

## Input parameters
MICROALGAE=$1
GENOME=$2
ANNOTATION=$3

## Building reference genome index
cd $GENOME
extract_splice_sites.py $ANNOTATION > splices_sites.ss
extract_exons.py $ANNOTATION > exons.exon
hisat2-build --ss splices_sites.ss --exon exons.exon $MICROALGAE.fa hisat2_index_$MICROALGAE

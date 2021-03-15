## Authors: Francisco J. Romero-Campero
##          Ana Belen Romero-Losada
## Contact: Francisco J. Romero-Campero - fran@us.es

#! /bin/bash

## Name input parameters
DATA=$1
PAIRED=$2
SAMPLE_FOLDER=$3
INDEX=$4
ANNOTATION=$5
NUM_SAMPLES=$6
INS=$7
EXP_FOLDER=$8
CONTROL=$9
EXPERIMENTAL=${10}
FOLD_CHANGE=${11}
Q_VALUE=${12}
FASTQ_LEFT=${13}
FASTQ_RIGHT=${14}
NPROC=${15}

ACC_NUMBER=${FASTQ_LEFT}

## Downloading or copying sample file depending on data source
cd ${SAMPLE_FOLDER} 
if [ $DATA == "DB" ]
then
	fastq-dump --split-files --gzip ${ACC_NUMBER}
elif [ $DATA == "FILES" ] && [ $PAIRED == "FALSE" ]
then
	cp ${ACC_NUMBER} sample_1.fastq.gz
	ACC_NUMBER=sample
elif [ $DATA == "FILES" ] && [ $PAIRED == "TRUE" ]
then
	cp ${FASTQ_LEFT} sample_1.fastq.gz
	cp ${FASTQ_RIGHT} sample_2.fastq.gz
	ACC_NUMBER=sample   
fi

## Sample quality control and read mapping to reference genome
if [ -f ${ACC_NUMBER}_2.fastq.gz ]
then
   #fastqc ${ACC_NUMBER}_1.fastq.gz
   #fastqc ${ACC_NUMBER}_2.fastq.gz

   hisat2 -p $NPROC --dta -x $INDEX -1 ${ACC_NUMBER}_1.fastq.gz -2 ${ACC_NUMBER}_2.fastq.gz -S ${ACC_NUMBER}.sam | tee mapping_stats
else
   fastqc ${ACC_NUMBER}_1.fastq.gz

   hisat2 -p $NPROC --dta -x $INDEX -U ${ACC_NUMBER}_1.fastq.gz -S ${ACC_NUMBER}.sam | tee mapping_stats
fi

## Generting sorted bam file
samtools sort -@ $NPROC -m 2G -o ${ACC_NUMBER}.bam ${ACC_NUMBER}.sam
rm ${ACC_NUMBER}.sam
rm *.fastq.gz
rm $HOME/ncbi/public/sra/${ACC_NUMBER}.sra
samtools index -@ $NPROC ${ACC_NUMBER}.bam
bamCoverage -p $NPROC -bs 10 --normalizeUsing CPM --bam ${ACC_NUMBER}.bam -o ${ACC_NUMBER}.bw
#rm *.bam

## Transcript assembly
stringtie -p $NPROC -G $ANNOTATION -o ${ACC_NUMBER}.gtf -l ${ACC_NUMBER} ${ACC_NUMBER}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLE_FOLDER}/${ACC_NUMBER}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -p $NPROC -e -B -G $ANNOTATION -o ${ACC_NUMBER}.gtf ${ACC_NUMBER}.bam

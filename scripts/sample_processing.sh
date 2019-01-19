## Author: Francisco J. Romero-Campero
## Date: January 2019
## Email: fran@us.es

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

## Name input parameters
SAMPLE_FOLDER=$1
ACC_NUMBER=$2
INDEX=$3
ANNOTATION=$4

## Downloading sample file
cd ${SAMPLE_FOLDER} 
#fastq-dump --split-files ${ACC_NUMBER}

## Sample quality control and read mapping to reference genome
#if [ -f ${ACC_NUMBER}_2.fastq ]
#then
#   fastqc ${ACC_NUMBER}_1.fastq
#   fastqc ${ACC_NUMBER}_2.fastq

#   hisat2 --dta -x $INDEX -1 ${ACC_NUMBER}_1.fastq -2 ${ACC_NUMBER}_2.fastq -S ${ACC_NUMBER}.sam
#else
#   fastqc ${ACC_NUMBER}_1.fastq

#   hisat2 --dta -x $INDEX -U ${ACC_NUMBER}_1.fastq -S ${ACC_NUMBER}.sam
#fi

## Generting sorted bam file
#samtools sort -o ${ACC_NUMBER}.bam ${ACC_NUMBER}.sam
#rm ${ACC_NUMBER}.sam
#rm *.fastq
#samtools index ${ACC_NUMBER}.bam
#bamCoverage -bs 10 --normalizeUsing CPM --bam ${ACC_NUMBER}.bam -o ${ACC_NUMBER}.bw

## Transcript assembly
stringtie -G $ANNOTATION -o ${ACC_NUMBER}.gtf -l ${ACC_NUMBER} ${ACC_NUMBER}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLE_FOLDER}/${ACC_NUMBER}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -e -B -G $ANNOTATION -o ${ACC_NUMBER}.gtf ${ACC_NUMBER}.bam


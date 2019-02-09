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
NUM_SAMPLES=$5
INS=$6
EXP_FOLDER=$7
CONTROL=$8
EXPERIMENTAL=$9
FOLD_CHANGE=${10}
Q_VALUE=${11}

## Downloading sample file
cd ${SAMPLE_FOLDER} 
fastq-dump --split-files ${ACC_NUMBER}

## Sample quality control and read mapping to reference genome
if [ -f ${ACC_NUMBER}_2.fastq ]
then
   fastqc ${ACC_NUMBER}_1.fastq
   fastqc ${ACC_NUMBER}_2.fastq

   hisat2 --dta -x $INDEX -1 ${ACC_NUMBER}_1.fastq -2 ${ACC_NUMBER}_2.fastq -S ${ACC_NUMBER}.sam
else
   fastqc ${ACC_NUMBER}_1.fastq

   hisat2 --dta -x $INDEX -U ${ACC_NUMBER}_1.fastq -S ${ACC_NUMBER}.sam
fi

## Generting sorted bam file
samtools -m 2G sort -o ${ACC_NUMBER}.bam ${ACC_NUMBER}.sam
rm ${ACC_NUMBER}.sam
rm *.fastq
rm $HOME/ncbi/public/sra/${ACC_NUMBER}.sra
samtools index ${ACC_NUMBER}.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam ${ACC_NUMBER}.bam -o ${ACC_NUMBER}.bw

## Transcript assembly
stringtie -G $ANNOTATION -o ${ACC_NUMBER}.gtf -l ${ACC_NUMBER} ${ACC_NUMBER}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLE_FOLDER}/${ACC_NUMBER}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -e -B -G $ANNOTATION -o ${ACC_NUMBER}.gtf ${ACC_NUMBER}.bam

## Write in blackboard
echo "SAMPLE " ${ACC_NUMBER} " DONE" >> ../../logs/blackboard.txt

## Count number of line in the blackboard to check the number of processed samples
PROCESSED_SAMPLES=$(wc -l ../../logs/blackboard.txt | awk '{print $1}')

## Submit scripts for transcriptome merging and differential gene expression 
if [ ${PROCESSED_SAMPLES} -eq ${NUM_SAMPLES} ]
then
   qsub -o ${EXP_FOLDER}/logs/transcriptome_merging $INS/scripts/transcriptome_merging.sh ${EXP_FOLDER} $ANNOTATION
   Rscript $INS/scripts/DE_analysis.R ${EXP_FOLDER}/samples $CONTROL $EXPERIMENTAL $FOLD_CHANGE $Q_VALUE
fi


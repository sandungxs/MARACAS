## Authors: Francisco J. Romero-Campero
##          Ana Belen Romero-Losada
## Contact: Francisco J. Romero-Campero - fran@us.es

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

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
EXPERIMENTAL=$10
FOLD_CHANGE=${11}
Q_VALUE=${12}
FASTQ_LEFT=${13}
FASTQ_RIGHT=${14}

ACC_NUMBER=${FASTQ_LEFT}

## Downloading or copying sample file depending on data source
cd ${SAMPLE_FOLDER} 
if [ $DATA == "DB" ]
then
	fastq-dump --split-files ${ACC_NUMBER}
elif [ $DATA == "FILES" ] && [ $PAIRED == "FALSE" ]
then
	cp ${ACC_NUMBER} sample_1.fastq.gz
	gunzip sample_1.fastq.gz
	ACC_NUMBER=sample
elif [ $DATA == "FILES" ] && [ $PAIRED == "TRUE" ]
then
	cp ${FASTQ_LEFT} sample_1.fastq.gz
	gunzip sample_1.fastq.gz
	cp ${FASTQ_RIGHT} sample_2.fastq.gz
	gunzip sample_2.fastq.gz
	ACC_NUMBER=sample   
fi

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
samtools sort -m 2G -o ${ACC_NUMBER}.bam ${ACC_NUMBER}.sam
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


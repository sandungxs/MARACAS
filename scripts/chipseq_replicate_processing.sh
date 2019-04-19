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
REPLICATE_FOLDER=$3
INDEX=$4
ANNOTATION=$5
NUM_REPLICATES=$6
INS=$7
CONTROL=$8
CHIP_LEFT=$9
CHIP_RIGHT=${10}
CONTROL_LEFT=${11}
CONTROL_RIGHT=${12}

## Downloading or copying sample file depending on data source
cd ${SAMPLE_FOLDER} 
if [ $DATA == "DB" ]
then

   fastq-dump --split-files ${CHIP_LEFT}
   if [ $PAIRED == "FALSE" ]
   then
      mv ${CHIP_LEFT}_1.fastq chip_1.fastq 
   elif [ $PAIRED == "TRUE" ]
   then
      mv ${CHIP_LEFT}_1.fastq chip_1.fastq
      mv ${CHIP_LEFT}_2.fastq chip_2.fastq
   fi

   if [ $CONTROL == "yes" ]
   then
      fastq-dump --split-files ${CONTROL_LEFT}
      if [ $PAIRED == "FALSE" ]
      then
         mv ${CONTROL_LEFT}_1.fastq control_1.fastq 
      elif [ $PAIRED == "TRUE" ]
      then
         mv ${CONTROL_LEFT}_1.fastq control_1.fastq
         mv ${CONTROL_LEFT}_2.fastq control_2.fastq
      fi
   fi

elif [ $DATA == "FILES" ] && [ $PAIRED == "FALSE" ]
then

   cp ${CHIP_LEFT} chip_1.fastq.gz
   gunzip chip_1.fastq.gz
   ACC_NUMBER=chip
   if [ $CONTROL == "yes" ]
   then
      cp ${CONTROL_LEFT} control_1.fastq.gz
      gunzip control_1.fastq.gz
   fi

elif [ $DATA == "FILES" ] && [ $PAIRED == "TRUE" ]
then

   cp ${CHIP_LEFT} chip_1.fastq.gz
   gunzip chip_1.fastq.gz
   cp ${CHIP_RIGHT} chip_2.fastq.gz
   gunzip chip_2.fastq.gz
   ACC_NUMBER=chip

   if [ $CONTROL == "yes" ]
   then
      cp ${CONTROL_LEFT} control_1.fastq.gz
      gunzip control_1.fastq.gz
      cp ${CONTROL_RIGHT} control_2.fastq.gz
      gunzip control_2.fastq.gz
   fi

fi

## Sample quality control, read mapping to reference genome, duplicates
## removal and bw generation. 
if [ $PAIRED == "yes" ]
then
   fastqc chip_1.fastq
   fastqc chip_2.fastq

   
   hisat2 --dta -x $INDEX -1 ${ACC_NUMBER}_1.fastq -2 ${ACC_NUMBER}_2.fastq -S ${ACC_NUMBER}.sam
else
   fastqc chip_1.fastq
   bowtie2 -x $INDEX -U chip_1.fastq -S raw_chip.sam
   samtools sort -n -m 2G -o raw_chip_sorted.bam chip.sam
   samtools fixmate -m raw_chip_sorted.bam chip_fixmate.bam
   samtools sort -m 2G -o chip_dup_sorted.bam chip_fixmate.bam 
   samtools markdup -r chip_dup_sorted.bam chip.bam
   samtools index chip.bam
   bamCoverage -bs 10 --normalizeUsing CPM --bam chip.bam -o chip.bw


   if [ $CONTROL == "yes" ]
   then
      bowtie2 -x $INDEX -U control_1.fastq -S raw_control.sam
      samtools sort -n -m 2G -o raw_chip_sorted.bam chip.sam
      samtools fixmate -m raw_chip_sorted.bam chip_fixmate.bam
      samtools sort -m 2G -o chip_dup_sorted.bam chip_fixmate.bam 
      samtools markdup -r chip_dup_sorted.bam chip.bam

   fi

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


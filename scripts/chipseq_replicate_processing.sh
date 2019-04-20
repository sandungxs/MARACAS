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
NUM_REPLICATES=$5
INS=$6
CONTROL=$7
CHIP_LEFT=$8
CHIP_RIGHT=${9}
CONTROL_LEFT=${10}
CONTROL_RIGHT=${11}

echo "DATA" $DATA
echo "PAIRED" $PAIRED
echo "replicate folder" ${REPLICATE_FOLDER}
echo "index" $INDEX
echo "number replicates" ${NUM_REPLICATES}
echo "instalation" $INS
echo "control" $CONTROL
echo "chip left" ${CHIP_LEFT}
echo "chip right" ${CHIP_RIGHT}
echo "control left" ${CONTROL_LEFT}
echo "control right" ${CONTROL_RIGHT}

## Downloading or copying sample file depending on data source
echo "access " ${REPLICATE_FOLDER}
echo "chip left " ${CHIP_LEFT}
echo "control left " ${CONTROL_LEFT}

cd ${REPLICATE_FOLDER}

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
   bowtie2 -x $INDEX -1 chip_1.fastq -2 chip_2.fastq -S raw_chip.sam
   if [ $CONTROL == "yes" ]
   then
      bowtie2 -x $INDEX -1 control_1.fastq -2 control_2.fastq -S raw_control.sam
   fi
else
   fastqc chip_1.fastq
   bowtie2 -x $INDEX -U chip_1.fastq -S raw_chip.sam

   if [ $CONTROL == "yes" ]
   then
      bowtie2 -x $INDEX -U control_1.fastq -S raw_control.sam
   fi
fi

## Duplicates removal and bw generation
samtools sort -n -m 2G -o raw_chip_sorted.bam raw_chip.sam
samtools fixmate -m raw_chip_sorted.bam chip_fixmate.bam
samtools sort -m 2G -o chip_dup_sorted.bam chip_fixmate.bam 
samtools markdup -r chip_dup_sorted.bam chip.bam
samtools index chip.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam chip.bam -o chip.bw

if [ $CONTROL == "yes" ]
then
   samtools sort -n -m 2G -o raw_control_sorted.bam raw_control.sam
   samtools fixmate -m raw_control_sorted.bam control_fixmate.bam
   samtools sort -m 2G -o control_dup_sorted.bam control_fixmate.bam 
   samtools markdup -r control_dup_sorted.bam control.bam
   samtools index control.bam
   bamCoverage -bs 10 --normalizeUsing CPM --bam control.bam -o control.bw
fi

## Remove files
rm *.sam
rm *.fastq
rm *fixmate*
rm *dup_sorted*
rm $HOME/ncbi/public/sra/${ACC_NUMBER}.sra

## Write in blackboard
echo "REPLICATE " ${ACC_NUMBER} " DONE" >> ../../logs/blackboard.txt

## Count number of line in the blackboard to check the number of processed samples
PROCESSED_REPLICATES=$(wc -l ../../logs/blackboard.txt | awk '{print $1}')

## Submit scripts for transcriptome merging and differential gene expression 
if [ ${PROCESSED_REPLICATES} -eq ${NUM_REPLICATES} ]
then
#   qsub -o ${EXP_FOLDER}/logs/transcriptome_merging $INS/scripts/transcriptome_merging.sh ${EXP_FOLDER} $ANNOTATION
#   Rscript $INS/scripts/DE_analysis.R ${EXP_FOLDER}/samples $CONTROL $EXPERIMENTAL $FOLD_CHANGE $Q_VALUE
   echo "Here Rscritp to be launched"
fi


#! /bin/bash

## Authors: Francisco J. Romero-Campero
##          Ana Belen Romero-Losada
## Contact: Francisco J. Romero-Campero - fran@us.es

## Name input parameters
DATA=$1
PAIRED=$2
WD=$3
MAIN_FOLDER=$4
INDEX=$5
NUM_REPLICATES=$6
INS=$7
CONTROL=$8
CHIP_LEFT=$9
CHIP_RIGHT=${10}
CONTROL_LEFT=${11}
CONTROL_RIGHT=${12}
CURRENT_REPLICATE=${13}
MODE=${14}
NPROC=${15}
ARCH=${16}
MICROALGAE=${17}
HISTONE=${18}
TF=${19}

REPLICATE_FOLDER=$WD/${MAIN_FOLDER}/replicates/replicate_${CURRENT_REPLICATE}

##echo "DATA" $DATA
echo "PAIRED" $PAIRED
##echo "replicate folder" ${REPLICATE_FOLDER}
##echo "index" $INDEX
##echo "number replicates" ${NUM_REPLICATES}
##echo "instalation" $INS
##echo "control" $CONTROL
##echo "chip left" ${CHIP_LEFT}
##echo "chip right" ${CHIP_RIGHT}
##echo "control left" ${CONTROL_LEFT}
##echo "control right" ${CONTROL_RIGHT}
echo "NPROC=" $NPROC

## Downloading or copying sample file depending on data source
##echo "access " ${REPLICATE_FOLDER}
##echo "chip left " ${CHIP_LEFT}
##echo "control left " ${CONTROL_LEFT}

cd ${REPLICATE_FOLDER}

if [ $DATA == "DB" ]
then

   echo ""
   echo "*******************************************"
   echo "* DOWNLOADING FASTQ FILES FOR CHIP SAMPLE *"
   echo "*******************************************"
   echo ""

   fastq-dump --gzip --split-files ${CHIP_LEFT}
   if [ $PAIRED == "FALSE" ]
   then
      mv ${CHIP_LEFT}_1.fastq.gz chip_1.fastq.gz 
   elif [ $PAIRED == "TRUE" ]
   then
      mv ${CHIP_LEFT}_1.fastq.gz chip_1.fastq.gz
      mv ${CHIP_LEFT}_2.fastq.gz chip_2.fastq.gz
   fi

   echo ""
   echo "********************************************"
   echo "* FASTQ FILES FOR CHIP SAMPLE DOWNLOADED!! *"
   echo "********************************************"
   echo ""

   if [ $CONTROL == "yes" ]
   then

      echo ""
      echo "********************************************"
      echo "* DOWNLOADING FASTQ FILES FOR INPUT SAMPLE *"
      echo "********************************************"
      echo ""

      fastq-dump --gzip --split-files ${CONTROL_LEFT}
      if [ $PAIRED == "FALSE" ]
      then
         mv ${CONTROL_LEFT}_1.fastq.gz control_1.fastq.gz 
      elif [ $PAIRED == "TRUE" ]
      then
         mv ${CONTROL_LEFT}_1.fastq.gz control_1.fastq.gz
         mv ${CONTROL_LEFT}_2.fastq.gz control_2.fastq.gz
      fi

      echo ""
      echo "**********************************************"
      echo "* FASTQ FILES FOR INPUT SAMPLE DOWNLOADED !! *"
      echo "**********************************************"
      echo ""

   fi

elif [ $DATA == "FILES" ] && [ $PAIRED == "FALSE" ]
then

   echo ""
   echo "***************************************"
   echo "* COPYING FASTQ FILES FOR CHIP SAMPLE *"
   echo "***************************************"
   echo ""

   cp ${CHIP_LEFT} chip_1.fastq.gz
   ACC_NUMBER=chip

   echo ""
   echo "**********************************************"
   echo "* COPYING FASTQ FILES FOR CHIP SAMPLE DONE!! *"
   echo "**********************************************"
   echo ""

   if [ $CONTROL == "yes" ]
   then

      echo ""
      echo "*****************************************"
      echo "* COPYING FASTQ FILES FOR INPUT SAMPLE  *"
      echo "*****************************************"
      echo ""

      cp ${CONTROL_LEFT} control_1.fastq.gz
      
      echo ""
      echo "************************************************"
      echo "* COPYING FASTQ FILES FOR INPUT SAMPLE DONE !! *"
      echo "************************************************"
      echo ""

   fi

elif [ $DATA == "FILES" ] && [ $PAIRED == "TRUE" ]
then

   echo ""
   echo "***************************************"
   echo "* COPYING FASTQ FILES FOR CHIP SAMPLE *"
   echo "***************************************"
   echo ""

   cp ${CHIP_LEFT} chip_1.fastq.gz
   cp ${CHIP_RIGHT} chip_2.fastq.gz
   ACC_NUMBER=chip

   echo ""
   echo "***********************************************"
   echo "* COPYING FASTQ FILES FOR CHIP SAMPLE DONE !! *"
   echo "***********************************************"
   echo ""

   if [ $CONTROL == "yes" ]
   then

      echo ""
      echo "*****************************************"
      echo "* COPYING FASTQ FILES FOR INPUT SAMPLE  *"
      echo "*****************************************"
      echo ""

      cp ${CONTROL_LEFT} control_1.fastq.gz
      cp ${CONTROL_RIGHT} control_2.fastq.gz
      
      echo ""
      echo "************************************************"
      echo "* COPYING FASTQ FILES FOR INPUT SAMPLE DONE !! *"
      echo "************************************************"
      echo ""

   fi
fi

echo ""
echo "********************************************************"
echo "* QUALITY CONTROL AND READ MAPPING TO REFERENCE GENOME *"
echo "********************************************************"
echo ""

## Sample quality control and read mapping to reference genome. 
if [ $PAIRED == "TRUE" ]
then
   fastqc chip_1.fastq.gz
   fastqc chip_2.fastq.gz
   bowtie2 -p $NPROC -x $INDEX -1 chip_1.fastq.gz -2 chip_2.fastq.gz -S raw_chip.sam &> chip_mapping_stats_${CURRENT_REPLICATE}
   if [ $CONTROL == "yes" ]
   then
      bowtie2 -p $NPROC -x $INDEX -1 control_1.fastq.gz -2 control_2.fastq.gz -S raw_control.sam &> control_mapping_stats_${CURRENT_REPLICATE}
   fi
else
   fastqc chip_1.fastq.gz
   bowtie2 -p $NPROC -x $INDEX -U chip_1.fastq.gz -S raw_chip.sam &> chip_mapping_stats_${CURRENT_REPLICATE}

   if [ $CONTROL == "yes" ]
   then
      fastqc control_1.fastq
      bowtie2 -p $NPROC -x $INDEX -U control_1.fastq.gz -S raw_control.sam &> control_mapping_stats_${CURRENT_REPLICATE}
   fi
fi

echo ""
echo "************************"
echo "* READ MAPPING DONE !! *"
echo "************************"
echo ""

echo ""
echo "************************************************"
echo "* PCR REMOVAL, BAM and BIGWIG FILES GENERATION *"
echo "************************************************"
echo ""

## Duplicates removal and bw generation
samtools sort -@ $NPROC -n -m 2G -o raw_chip_sorted.bam raw_chip.sam
samtools fixmate -@ $NPROC -m raw_chip_sorted.bam chip_fixmate.bam
samtools sort -@ $NPROC -m 2G -o chip_dup_sorted.bam chip_fixmate.bam 
samtools markdup -@ $NPROC -r chip_dup_sorted.bam chip.bam
samtools index -@ $NPROC chip.bam
bamCoverage -p $NPROC -bs 5 --normalizeUsing CPM --bam chip.bam -o chip_${CURRENT_REPLICATE}.bw

if [ $CONTROL == "yes" ]
then
   samtools sort -@ $NPROC -n -m 2G -o raw_control_sorted.bam raw_control.sam
   samtools fixmate -@ $NPROC -m raw_control_sorted.bam control_fixmate.bam
   samtools sort -@ $NPROC -m 2G -o control_dup_sorted.bam control_fixmate.bam 
   samtools markdup -@ $NPROC -r control_dup_sorted.bam control.bam
   samtools index -@ $NPROC control.bam
   bamCoverage -p $NPROC -bs 5 --normalizeUsing CPM --bam control.bam -o control_${CURRENT_REPLICATE}.bw
fi

## Remove files
rm *.sam
rm *.fastq*
rm *fixmate*
rm *dup_sorted*
rm raw*

echo ""
echo "*******************************"
echo "* BIGWIG FILES GENERATED !!! *"
echo "******************************"
echo ""

echo ""
echo "****************"
echo "* PEAK CALLING *"
echo "****************"
echo ""

## Peak calling 
if [ $MODE == "histone_modification" ]
then

   if [ $CONTROL == "yes" ]
   then
      echo "Peak calling with control sample"
      macs2 callpeak -t chip.bam -c control.bam -f BAM --outdir . -n replicate_${CURRENT_REPLICATE} --nomodel &> macs_output
   else
      echo "Peak calling without control sample"
      macs2 callpeak -t chip.bam -f BAM --outdir. -n replicate_${CURRENT_REPLICATE} --nomodel  &> macs_output
   fi

elif [ $MODE == "transcription_factor" ]
then
   if [ $CONTROL == "yes" ]
   then
      echo "Peak calling with control sample"
      macs2 callpeak -t chip.bam -c control.bam -f BAM --outdir . -n replicate_${CURRENT_REPLICATE} &> macs_output
   else
      echo "Peak calling without control sample"
      macs2 callpeak -t chip.bam -f BAM --outdir. -n replicate_${CURRENT_REPLICATE} &> macs_output
   fi

else
   echo "Incorrect value for MODE " $MODE
   echo "Only histone_modification or transcription_factor are allowed"
   exit
fi

rm *.bam*

echo ""
echo "************************"
echo "* PEAK CALLING DONE !!!*"
echo "************************"
echo ""

if [ $ARCH == "SLURM" ]
then
   ## Write in blackboard
   echo "REPLICATE " ${CURRENT_REPLICATE} " DONE" >> $WD/${MAIN_FOLDER}/logs/blackboard.txt

   ## Count number of line in the blackboard to check the number of processed samples
   PROCESSED_REPLICATES=$(wc -l $WD/${MAIN_FOLDER}/logs/blackboard.txt | awk '{print $1}')

   ## Submit scripts for transcriptome merging and differential gene expression 
   if [ ${PROCESSED_REPLICATES} -eq ${NUM_REPLICATES} ]
   then
      if [ $MODE == "histone_modification" ]
      then
         sbatch $INS/scripts/chipseq_final_processing.sh $WD ${MAIN_FOLDER} ${NUM_REPLICATES} \
                                                          $MODE ${INS} $MICROALGAE $HISTONE $HISTONE \
                                                          $CONTROL $NPROC
      elif [ $MODE == "transcription_factor" ]
      then
         sbatch $INS/scripts/chipseq_final_processing.sh $WD ${MAIN_FOLDER} ${NUM_REPLICATES} \
                                                          $MODE ${INS} $MICROALGAE $TF $TF \
                                                          $CONTROL $NPROC
      fi
   fi
fi

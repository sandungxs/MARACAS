#! /bin/bash

## Authors: Francisco J. Romero-Campero
##          Ana Belen Romero-Losada
## Contact: Francisco J. Romero-Campero - fran@us.es

## Name input parameters
DATA=$1
PAIRED=$2
SAMPLE_FOLDER=$3
INDEX=$4
ANNOTATION=$5
NUM_SAMPLES=$6
EXP_FOLDER=$7
CONTROL=$8
EXPERIMENTAL=$9
FOLD_CHANGE=${10}
Q_VALUE=${11}
FASTQ_LEFT=${12}
FASTQ_RIGHT=${13}
NPROC=${14}
ARCH=${15}
MICROALGAE=${16}
MAPPER=${17}

ACC_NUMBER=${FASTQ_LEFT}

## Downloading or copying sample file depending on data source
echo ""
echo "* Downloading or copying files *" 
echo "*******************************"
echo ""

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
	ACC_NUMBER=sample   ####y luego como diferencio una muestra de otra? 
fi

## Sample quality control and read mapping to reference genome
if [ $MAPPER == "hisat2" ]
then
if [ -f ${ACC_NUMBER}_2.fastq.gz ]
then

   echo ""
   echo "* Quality Analysis *" 
   echo "********************"
   echo ""
   
   fastqc ${ACC_NUMBER}_1.fastq.gz
   fastqc ${ACC_NUMBER}_2.fastq.gz
   
   echo ""
   echo "*    Read Mapping  *" 
   echo "********************"
   echo ""

   hisat2 -p $NPROC --dta -x $INDEX -1 ${ACC_NUMBER}_1.fastq.gz -2 ${ACC_NUMBER}_2.fastq.gz -S ${ACC_NUMBER}.sam --summary-file mapping_stats
else

   echo ""
   echo "* Quality Analysis *" 
   echo "********************"
   echo ""

   fastqc ${ACC_NUMBER}_1.fastq.gz
   
   echo ""
   echo "* Read Mapping     *" 
   echo "********************"
   echo ""
   
   hisat2 -p $NPROC --dta -x $INDEX -U ${ACC_NUMBER}_1.fastq.gz -S ${ACC_NUMBER}.sam --summary-file mapping_stats
fi

echo ""
echo "* Generating mapping signal in bigwig format *" 
echo "**********************************************"
echo ""

## Generting sorted bam file
samtools sort -@ $NPROC -m 2G -o ${ACC_NUMBER}.bam ${ACC_NUMBER}.sam
rm ${ACC_NUMBER}.sam
#rm *.fastq.gz
#rm $HOME/ncbi/public/sra/${ACC_NUMBER}.sra
samtools index -@ $NPROC ${ACC_NUMBER}.bam
bamCoverage -p $NPROC -bs 10 --normalizeUsing CPM --bam ${ACC_NUMBER}.bam -o ${ACC_NUMBER}.bw

## Transcript assembly
echo ""
echo "* Transcript assembly and gene expression estimation *" 
echo "******************************************************"
echo ""

stringtie -p $NRPOC -G $ANNOTATION -o ${ACC_NUMBER}.gtf -l ${ACC_NUMBER} ${ACC_NUMBER}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLE_FOLDER}/${ACC_NUMBER}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -p $NPROC -e -B -G $ANNOTATION -o ${ACC_NUMBER}.gtf ${ACC_NUMBER}.bam
rm ${ACC_NUMBER}.bam

sed -i 's/#/HHAASSHH/g' t_data.ctab

## Synchronization
if [ $ARCH == "SLURM" ]
then
   ## Write in blackboard
   echo "SAMPLE " ${CURRENT_REPLICATE} " DONE" >> ${SAMPLE_FOLDER}/../../logs/blackboard.txt

   ## Count number of line in the blackboard to check the number of processed samples
   PROCESSED_SAMPLES=$(wc -l ${SAMPLE_FOLDER}/../../logs/blackboard.txt | awk '{print $1}')

   ## Submit scripts for transcriptome merging and differential gene expression 
   if [ ${PROCESSED_SAMPLES} -eq ${NUM_SAMPLES} ]
   then
      
      sbatch $MARACAS/scripts/transcriptome_merging.sh ${SAMPLE_FOLDER}/../../ $MARACAS/data/${MICROALGAE}/annotation/${MICROALGAE}.gtf
      
      echo ""
      echo "* Computing differential gene expression *" 
      echo "******************************************"
      echo ""
      Rscript $MARACAS/scripts/DE_analysis.R ${SAMPLE_FOLDER}/../ ${CONTROL} ${EXPERIMENTAL} $FOLD_CHANGE $Q_VALUE $MICROALGAE
      
      echo ""
      echo "* Generating output reports *" 
      echo "*****************************"
      echo ""

      Rscript -e "rmarkdown::render('${SAMPLE_FOLDER}/../../results/DE_report.Rmd', 'pdf_document')" 
      Rscript -e "rmarkdown::render('${SAMPLE_FOLDER}/../../results/DE_report.Rmd', 'html_document')" 

      sed -i 's/HHAASSHH/#/g' ${SAMPLE_FOLDER}/../../results/*
   fi
fi

fi


############################# Sample processing with kallisto##################

if [ $MAPPER == "kallisto" ]
then

## Downloading or copying sample file depending on data source
echo ""
echo "* Downloading or copying files *" 
echo "*******************************"
echo ""

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

if [ -f ${ACC_NUMBER}_2.fastq.gz ]
then

   echo ""
   echo "* Quality Analysis *" 
   echo "********************"
   echo ""
   
   fastqc ${ACC_NUMBER}_1.fastq.gz
   fastqc ${ACC_NUMBER}_2.fastq.gz
   
   echo ""
   echo "*    Read Mapping  *" 
   echo "********************"
   echo ""

    kallisto quant -i ${MAPPER}_index_${MICROALGAE}.idx  -b 100 -o kallisto_out_${ACC_NUMBER} --genomebam --plaintext --gtf ${ANNOTATION} --chromosomes $MARACAS/data/${MICROALGAE}/genome/chrom_${MICROALGAE} -l 200 -s 10 --fr-stranded ${ACC_NUMBER}_1.fastq.gz --rf-stranded ${ACC_NUMBER}_2.fastq.gz

else

   echo ""
   echo "* Quality Analysis *" 
   echo "********************"
   echo ""

   fastqc ${ACC_NUMBER}_1.fastq.gz
   
   echo ""
   echo "* Read Mapping     *" 
   echo "********************"
   echo ""
   
    kallisto quant -i ${MAPPER}_index_${MICROALGAE}.idx  -b 100 -o kallisto_out_${ACC_NUMBER} --genomebam --plaintext --gtf ${ANNOTATION} --chromosomes $MARACAS/data/${MICROALGAE}/genome/chrom_${MICROALGAE} -l 200 -s 10 --single ${ACC_NUMBER}_1.fastq.gz 
fi

## Synchronization
if [ $ARCH == "SLURM" ]
then
   ## Write in blackboard
   echo "SAMPLE " ${CURRENT_REPLICATE} " DONE" >> ${SAMPLE_FOLDER}/../../logs/blackboard.txt

   ## Count number of line in the blackboard to check the number of processed samples
   PROCESSED_SAMPLES=$(wc -l ${SAMPLE_FOLDER}/../../logs/blackboard.txt | awk '{print $1}')

   ## Submit scripts for transcriptome merging and differential gene expression 
   if [ ${PROCESSED_SAMPLES} -eq ${NUM_SAMPLES} ]
   then
      
      sbatch $MARACAS/scripts/transcriptome_merging.sh ${SAMPLE_FOLDER}/../../ $MARACAS/data/${MICROALGAE}/annotation/${MICROALGAE}.gtf
      
      echo ""
      echo "* Computing differential gene expression *" 
      echo "******************************************"
      echo ""
      Rscript $MARACAS/scripts/DE_analysis.R ${SAMPLE_FOLDER}/../ ${CONTROL} ${EXPERIMENTAL} $FOLD_CHANGE $Q_VALUE $MICROALGAE $MAPPER
      
      echo ""
      echo "* Generating output reports *" 
      echo "*****************************"
      echo ""

      Rscript -e "rmarkdown::render('${SAMPLE_FOLDER}/../../results/DE_report.Rmd', 'pdf_document')" 
      Rscript -e "rmarkdown::render('${SAMPLE_FOLDER}/../../results/DE_report.Rmd', 'html_document')" 

      sed -i 's/HHAASSHH/#/g' ${SAMPLE_FOLDER}/../../results/*
   fi
fi

fi #final

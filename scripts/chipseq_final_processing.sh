#! /bin/bash 

WD=$1
MAIN_FOLDER=$2
NUM_REPLICATES=$3
MODE=$4
MICROALGAE=$5 
HISTONE=$6 
TF=$7
INCLUDED_CONTROL=$8
NPROC=$9

echo "------------------------"
echo "WD=" $WD
echo "MAIN_FOLDER=" ${MAIN_FOLDER}
echo "NUM_REPLICATES=" ${NUM_REPLICATES}
echo "MODE=" $MODE
echo "MICRO=" $MICROALGAE
echo "HISTONE=" $HISTONE
echo "TF=" $TF
echo "INCLUDED_CONTROL=" ${INCLUDED_CONTROL}
echo "NPROC=" $NPROC
echo "----------------------"

echo ""
echo "---------------------------------------------"
echo "---------------------------------------------"
echo "|| STEP 4: FINAL RESULT REPORT GENERATION. ||"
echo "---------------------------------------------"
echo "---------------------------------------------"
echo ""

if [ ${NUM_REPLICATES} -gt 1 ]
then
   cd $WD/${MAIN_FOLDER}/results
   intersectBed -a  $WD/${MAIN_FOLDER}/replicates/replicate_1/replicate_1_peaks.narrowPeak \
                -b  $WD/${MAIN_FOLDER}/replicates/replicate_2/replicate_2_peaks.narrowPeak > acum_peaks.narrowPeak
   bigwigCompare -p $NPROC -b1 $WD/${MAIN_FOLDER}/replicates/replicate_1/chip_1.bw \
                 -b2 $WD/${MAIN_FOLDER}/replicates/replicate_2/chip_2.bw \
                 --scaleFactors 1:1 --operation add \
                 --binSize 5 --outFileName sum_chip.bw \
                 --outFileFormat bigwig 
   for  i in `seq 3 ${NUM_REPLICATES}`
   do
      intersectBed -a acum_peaks.narrowPeak -b ../replicates/replicate_$i/replicate_${i}_peaks.narrowPeak > acum_peaks_2.narrowPeak
      rm acum_peaks.narrowPeak
      mv acum_peaks_2.narrowPeak acum_peaks.narrowPeak

      cp sum_chip.bw cp_sum_chip.bw
      rm sum_chip.bw
      bigwigCompare -p $NPROC -b1 ../replicates/replicate_$i/chip_$i.bw -b2 cp_sum_chip.bw \
                    --scaleFactors 1:1 --operation add \
                    --binSize 5 --outFileName sum_chip.bw \
                    --outFileFormat bigwig
      rm cp_sum_chip.bw 

   done
   mv acum_peaks.narrowPeak output_peaks.narrowPeak

   ## Add bc dependency
   SCALEFACTOR=$(echo "1/(2*${NUM_REPLICATES})" | bc -l)
   bigwigCompare -p $NPROC -b1 sum_chip.bw -b2 sum_chip.bw \
                 --scaleFactors $SCALEFACTOR:$SCALEFACTOR --operation add \
                 --binSize 5 --outFileName chip_average.bw \
                 --outFileFormat bigwig 
   
else
   cp $WD/${MAIN_FOLDER}/replicates/replicate_1_peaks.narrowPeak $WD/${MAIN_FOLDER}/results/output_peaks.narrowPeak
fi

if [ $MODE == "histone_modification" ]
then
   Rscript $MARACAS/scripts/create_Rmd.R $WD/${MAIN_FOLDER}/results/ChIP_seq_analysis_report.Rmd $MICROALGAE $MODE $HISTONE ${NUM_REPLICATES} ${INCLUDED_CONTROL}
elif [ $MODE == "transcription_factor" ]
then
   Rscript $MARACAS/scripts/create_Rmd.R $WD/${MAIN_FOLDER}/results/ChIP_seq_analysis_report.Rmd $MICROALGAE $MODE $TF ${NUM_REPLICATES} ${INCLUDED_CONTROL}
fi

# Add rmarkdown dependency
Rscript -e "rmarkdown::render('ChIP_seq_analysis_report.Rmd', 'pdf_document')"
Rscript -e "rmarkdown::render('ChIP_seq_analysis_report.Rmd', 'html_document')"

echo ""
echo ""
echo "**************************************************************************************************************"
echo "*                                    ANALYSIS DONE!!!                                                        *"
echo "* OUTPUT REPORTS IN HTML AND PDF FORMATS HAVE BEEN GENERATED IN THE FOLDER $WD/${MAIN_FOLDER}/results"
echo "*                                  ENJOY YOUR RESULTS!!!                                                     *"
echo "**************************************************************************************************************"
echo ""
echo ""
echo ""
echo "                                    ░░░░░░░░░░░░░░░░░░░░░░█████████"
echo "                                    ░░███████░░░░░░░░░░███▒▒▒▒▒▒▒▒███"
echo "                                    ░░█▒▒▒▒▒▒█░░░░░░░███▒▒▒▒▒▒▒▒▒▒▒▒▒███"
echo "                                    ░░░█▒▒▒▒▒▒█░░░░██▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒██"
echo "                                    ░░░░█▒▒▒▒▒█░░░██▒▒▒▒▒██▒▒▒▒▒▒██▒▒▒▒▒███"
echo "                                    ░░░░░█▒▒▒█░░░█▒▒▒▒▒▒████▒▒▒▒████▒▒▒▒▒▒██"
echo "                                    ░░░█████████████▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒██"
echo "                                    ░░░█▒▒▒▒▒▒▒▒▒▒▒▒█▒▒▒▒▒▒▒▒▒█▒▒▒▒▒▒▒▒▒▒▒██"
echo "                                    ░██▒▒▒▒▒▒▒▒▒▒▒▒▒█▒▒▒██▒▒▒▒▒▒▒▒▒▒██▒▒▒▒██"
echo "                                    ██▒▒▒███████████▒▒▒▒▒██▒▒▒▒▒▒▒▒██▒▒▒▒▒██"
echo "                                    █▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒█▒▒▒▒▒▒████████▒▒▒▒▒▒▒██"
echo "                                    ██▒▒▒▒▒▒▒▒▒▒▒▒▒▒█▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒██"
echo "                                    ░█▒▒▒███████████▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒██"
echo "                                    ░██▒▒▒▒▒▒▒▒▒▒████▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒█"
echo "                                    ░░████████████░░░█████████████████"
echo ""
echo ""
echo ""





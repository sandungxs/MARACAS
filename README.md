# MARACAS - MicroAlgae RnA-seq and Chip-seq AnalysiS



[![DOI](https://zenodo.org/badge/165233815.svg)](https://zenodo.org/badge/latestdoi/165233815) [![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


## What is MARACAS?

MARACAS is an automatic computational pipeline specifically designed for the analysis of **microalgae** RNA-seq and ChIP-seq
data. MARACAS starts processing raw fastq files and generates lists of differentially 
expressed genes for RNA-seq data and lists of genomic loci in bed format for ChIP-seq data. BigWig files containing normalized mapping signal are also generated. The analysis are performed according to user specified parameters. 
Reports in html and pdf format are produced for easy exploration of the results. 

Differential expressed gene lists and genomic loci lists can be further analyzed using [**our online tool AlgaeFUN for
functional analysis**](https://greennetwork.us.es/AlgaeFUN/).

MARACAS supports a wide range of microalge including:

* Ostreococcus tauri
* Micromonas pusilla CCMP1545
* Bathycoccus prasinos
* Coccomyxa subellipsoidea
* Chlamydomonas reinhardtii
* Volvox Carteri
* Dunaliella salina
* Haematococcus lacustris
* Chlomochloris zofingiensis
* Klebsormidium nitens
* Mesotaenium endlicherianum
* Spirogloea muscicola
* Phaeodactylum tricornutum
* Nannochloropsis gaditana

MARACAS can be executed in a sequential mode in a laptop or server and in a distributed/parallel mode in a computer cluster.

## How to install it?

MARACAS requires the following dependencies that need to be previously installed:

* [FASTQC: A quality control tool for high throughput sequence data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [HISAT2: graph-based alignment of next generation sequencing reads to a population of genomes](http://daehwankimlab.github.io/hisat2/)
* [Bowtie 2: an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format](http://www.htslib.org/)
* [deepTools: tools for exploring deep sequencing data](https://deeptools.readthedocs.io/en/develop/index.html)
* [StringTie: Transcript assembly and quantification for RNA-Seq](https://ccb.jhu.edu/software/stringtie/)
* [MACS2: Model-based Analysis for ChIP-Seq](https://pypi.org/project/MACS2/)
* [R: free software environment for statistical computing and graphics with the following packages:](https://www.r-project.org/)
    * [ballgown: Tools for statistical analysis of assembled transcriptomes, including flexible differential expression analysis, visualization of transcript structures, and matching of assembled transcripts to annotation.](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html)
    * [FactoMineR: an R package dedicated to multivariate Exploratory Data Analysis.](http://factominer.free.fr/)
    * [factoextra: an R package making easy to extract and visualize the output of exploratory multivariate data analyses](https://cran.r-project.org/web/packages/factoextra/index.html)
    * [limma: Linear Models for Microarray and RNA-seq Data](https://bioconductor.org/packages/release/bioc/html/limma.html)
    * [rmarkdown: Dynamic Documents for R, Convert R Markdown documents into a variety of formats.](https://cran.r-project.org/web/packages/rmarkdown/index.html)
* [TeX Live: TeX document production system](https://www.tug.org/texlive/)


To install MARACAS you do NOT need sudo/administrator permissions just follow these steps:

1. Download the code from Github, for example, to your opt folder (make sure you have such a folder or create it with mkdir opt in your home directory):

```
cd
mkdir opt
cd opt
git clone https://github.com/fran-romero-campero/MARACAS.git
```

2. Add to your PATH variable defined in your .bashrc the path to the scripts folder so the MARACAS scripts can be executed from the command line: 

```
cd
echo "PATH=$PATH:$HOME/opt/MARACAS/scripts" >> .bashrc
source .bashrc
```

3. Add to your .bashrc a new variable MARACAS indicating the path to the MARACAS folder:  

```
echo "export MARACAS=$HOME/opt/MARACAS/" >> .bashrc
source .bashrc
```

## How to use it?

To run our pipeline MARACAS you first have to create a parameter file. This parameter file is different for RNA-seq or ChIP-seq data analysis. You can find examples of these files in the **test** and **examples** folders. 

### Run MARACAS and Parameter file for RNA-seq data Analysis

To process RNA-seq data use the executable **maracas-rna-seq** with a single input consisting in the parameter file:

```
maracas-rna-seq <parameter-file-rna-seq>
```

[Click here for an example of the parameter file used for RNA-seq analysis in sequential mode.](https://github.com/fran-romero-campero/MARACAS/blob/master/test/rna_seq_test_params_file_sequential.txt)

[Click here for an example of the parameter file used for RNA-seq analysis in distributed mode using SLURM.](https://github.com/fran-romero-campero/MARACAS/blob/master/test/rna_seq_test_params_file_distributed.txt)

Next the parameters to be specified in the parameter file for RNA-seq analysis are listed:

* **data_source:** This parameter specifies the source of the data to be analysed. It can take the value **FILES** when the fastq files are already located in a folder in the computer where MARACAS is installed or the value **DB** when the data has been freely deposited in the GEO data base.

* **cluster:** This parameter specifies the execution mode for MARACAS. It can take two different values: **SERVER** or **SLURM**. When SERVER mode is specified MARACAS will be 
executed with a sequential analysis of the different samples. Whereas when SLURM is specified MARACAS will be executed in a parallel/distributed manner processing samples simultaneously in different computational nodes. In this last case SLURM needs to be installed in your computer cluster. 

* **number_processors:** This parameter specifies the number of processors that can be
used by MARACAS. It can take an integer value for instance 2.

* **paired_end:** This parameter can take the values **FALSE** when your data is single end (a single fastq file per sample) and **TRUE** when your data is paired end (two fastq files R1 and R2 per sample).

* **working_directory:** This parameter specifies the directory where the output folder 
containing the results of the RNA-seq data analysis will be generated. For example, */home/cool_user/research/*

* **microalgae:** This parameter specifies the name of the microalge from which the
RNA-seq data was generated. It can take one the following exact names:
 
    * ostreococcus_tauri
    * micromonas_pusillaCCMP1545
    * bathycoccus_prasinos
    * coccomyxa_subellipsoidea
    * chlamydomonas_reinhardtii
    * volvox_carteri
    * dunaliella_salina
    * haematococcus_lacustris
    * chlomochloris_zofingiensis
    * klebsormidium_nitens
    * mesotaenium_endlicherianum
    * spirogloea_muscicola
    * phaeodactylum_tricornutum
    * nannochloropsis_gaditana


* **main_folder:** This parameter specifies the name of the folder where the 
final and intermediary results of the analysis will be saved. 

* **number_of_samples:** This parameter specifies the total number of samples to be analyzed. 

* **control_condition_name:** This parameter specifies the name of the control condition.

* **experimental_condition_name:** This parameter specifies the name of the experimental condition.

* **loc_sampleN:** In the case of single end data (*paired_end: FALSE*) and when the data to be processed are stored in your computer in fastq format (*data_source: FILES*) this parameter specifies the path and file name of sampleN where N=1,...,number_of_samples.

* **acc_sampleN:** In the case when the data to be processed needs to be retrieve from the GEO database (*data_source: DB*) this parameter specifies accession number identifying the fastq file in GEO.

* **loc_sample_leftN:** In the case of paired end data (*paired_end: TRUE*) and when the data to be processed are stored in your computer in fastq format (*data_source: FILES*) this parameter specifies the path and name of the fastq file containing the left reads for sampleN where N=1,...,number_of_samples.

* **loc_sample_rightN:** In the case of paired end data (*paired_end: TRUE*) and when the data to be processed are stored in your computer in fastq format (*data_source: FILES*) this parameter specifies the path and name of the fastq file containing the right reads for sampleN where N=1,...,number_of_samples.

* **condition_sample1:** This parameter specifies for each sampleN where N=1,...,number_of_samples the condition name to which it corresponds. The condition name has to be one of the names specified in the previous parameters *control_condition_name* or *experimental_condition_name*. 

* **fold_change:** This parameter specifies the fold-change used to determine differential expressed genes in the experimental condition when compared to the control condition. 

* **q_value:** This parameter specifies the q-value, adjusted p-value or FDR used to determine differential expressed genes in the experimental condition when compared to the control condition. 


### Run MARACAS and Parameter file for ChIP-seq data Analysis 

To process ChIP-seq data use the executable **maracas-chip-seq** with a single input consisting in the parameter file:

```
maracas-chip-seq <parameter-file-rna-seq>
```
[Click here for an example of the parameter file used for ChIP-seq analysis in sequential mode.](https://github.com/fran-romero-campero/MARACAS/blob/master/test/chip_seq_test_params_file_sequential.txt)

[Click here for an example of the parameter file used for ChIP-seq analysis in distributed mode using SLURM.](https://github.com/fran-romero-campero/MARACAS/blob/master/test/chip_seq_test_params_file_distributed.txt)

The parameters **data_source:**, **cluster:**, **number_processors:**, **working_directory:**, **microalgae:**, **main_folder:**, **number_of_replicates:** and **paired_end:** are the same as for the case of an RNA-seq analysis. You can find their description in the previous section.

Next the specific parameters for a ChIP-seq analysis are listed:

* **included_control:** This parameter can take the values *yes* or *no* depending whether or not your experimental design includes a control condition such as an input, mock or similar. 

* **mode:** This parameter can take the values *transcription_factor* or *histone_modification* to specify if your ChIP-seq data was generated for a transcription factor or a histone modification. 

* **transcription_factor:** When *mode: transcription_factor* this parameter specifies the name of the corresponding transcription factor.

* **histone_modification:** When *mode: histone_modification* this parameter specifies the name of the corresponding histone modification.

* **chip_replicate_N:** When *data_source: DB* this parameter specifies the accession number in GEO corresponding to replicate N of the ChIP sample.   

* **control_replicate_N:** When *data_source: DB* and *included_control: yes* this parameter specifies the accession number in GEO corresponding to replicate N of the control sample.

* **loc_chip_replicate_N:** When *data_source: FILES* and *paired_end: FALSE* this parameter specifies the path including the fastq file name corresponding to replicate N of the ChIP sample.  

* **loc_control_replicate_N:** When *data_source: FILES*, *included_control: yes* and *paired_end: FALSE* this parameter specifies the path including the file name corresponding to replicate N of the control sample.

* **loc_chip_replicate_left_N:**, **loc_chip_replicate_right_N:**, **loc_control_replicate_left_N:** and **loc_control_replicate_right_N:** These parameters specify the path including the fastq file names corresponding to the chip and control samples when *data_source: FILES*, *included_control: yes* and *paired_end: TRUE*.    

### Results generated by MARACAS 

The final reports and the intermediary result files generated during the execution of the MARACAS pipeline will be saved in the folder specified in the parameter *main_folder*. 

* [Click here to download a compressed example folder generated by MARACAS for an RNA-seq analysis.](https://github.com/fran-romero-campero/MARACAS/blob/master/examples/rna_seq_test.zip)
* [Click here to download a compressed example folder generated by MARACAS for a ChIP-seq analysis.](https://github.com/fran-romero-campero/MARACAS/blob/master/examples/chip_test.zip)

The MARACAS output folder for an RNA-seq data analysis consists of the following subfolders:

* **samples** This folder contains a subfolder for each sample whith files for quality control, mapping stats, normalized mapping signal in bigwig format, transcript assembly and mapping information for exons, transcripts and introns. 

* **results** This folder contains the information of the final results collected in two reports in html and pdf format. Click here for an example of result report in pdf format and an example of report in html format.

## Test Data Sets

In order to test that MARACAS is correctly installed in your computer or server make sure you have a **tmp** folder in your home directory and run the following instructions:

* A test for RNA-seq data analysis in your computer or server:

```
maracas-rna-seq $MARACAS/test/rna_seq_test_params_file_sequential.txt
```

* A test for ChIP-seq data analysis in your computer or server:

```
maracas-chip-seq $MARACAS/test/chip_seq_test_params_file_sequential.txt
```

In order to test that MARACAS is correctly installed in your cluster with SLURM make sure you have a **tmp** folder in your home directory and run the following instructions:

* A test for RNA-seq data analysis in your SLURM cluster:

```
maracas-rna-seq $MARACAS/test/rna_seq_test_params_file_distributed.txt
```
* A test for ChIP-seq data analysis in your SLURM cluster:

```
maracas-chip-seq $MARACAS/test/chip_seq_test_params_file_distributed.txt
```


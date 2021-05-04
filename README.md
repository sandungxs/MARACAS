# MARACAS - MicroAlgae RnA-seq and Chip-seq AnalysiS

## What is MARACAS?

MARACAS is an automatic computational pipeline specifically designed for the analysis of RNA-seq and ChIP-seq
data for **microalgae**. MARACAS starts processing raw fastq files and it generates lists of differentially 
expressed genes for RNA-seq data and lists of genomic loci in bed format for ChIP-seq data. BigWig files with the 
normalized mapping signal are also generated. The analysis are performed according to user specified parameters. 
Reports in html and pdf format are produced for easy exploration of the results. 

Differential expressed gene lists and genomic loci lists can be further analyzed using our online tool AlgaeFUN for
functional analysis.

MARACAS supports a wide range of microalge including:

* Ostreococcus tauri
* Chlamydomonas reinhardtii
* Haematococcus lacustris
* Dunaliella salina
* Volvox Carteri
* Phaeodactylum tricornutum
* Nannochloropsis gaditana
* Ostreococcus lucimarinus
* Coccomyxa subellipsoidea
* Bathycoccus prasinos
* Klebsormidium nitens
* Chlomochloris zofingiensis
* Micromonas pusilla CCMP1545

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


To install MARACAS you do NOT need sudo/operator/administrator permissions just follow these steps:

1. Download the code from Github and uncompress it, for example, to your opt folder:

```
cd opt
wget link
unzip master
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

To run our pipeline MARACAS you first have to create a parameter file. This parameter file is different for RNA-seq or ChIP-seq data analysis. You can find examples of these files in the **test** folder. In the **test** folder you can also find a template for each type of parameter file. 

### Run MARACAS and Parameter file for RNA-seq data Analysis

To process RNA-seq data use the executable **maracas-rna-seq** with a single input consisting in the parameter file:

```
maracas-rna-seq <parameter-file-rna-seq>
```

Click here for an example of the parameter file used for RNA-seq analysis. 
Click here for a template of the parameter file used for RNA-seq analysis. 

Next the parameters to be specified in the parameter file for RNA-seq analysis is listed:

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

    * bathycoccus_prasinos
    * chlamydomonas_reinhardtii
    * chromochloris_zofingiensis
    * coccomyxa_subellipsoidea
    * dunaliella_salina
    * klebsormidium_nitens
    * micromonas_pusillaCCMP1545
    * nannochloropsis_gaditana
    * ostreococcus_tauri
    * phaeodactylum_tricornutum
    * raphidocelis_subcapitata
    * volvox_carteri

* **main_folder:** This parameter specifies the name of the folder where the 
output local_test
number_of_samples: 4
control_condition_name: iron
experimental_condition_name: no_iron
loc_sample1: /home/fran/research/MARACAS/test/sample_1.fastq.gz
condition_sample1: iron
loc_sample2: /home/fran/research/MARACAS/test/sample_2.fastq.gz
condition_sample2: iron
loc_sample3: /home/fran/research/MARACAS/test/sample_3.fastq.gz
condition_sample3: no_iron
loc_sample4: /home/fran/research/MARACAS/test/sample_4.fastq.gz
condition_sample4: no_iron
fold_change: 2
q_value: 0.5


### Run MARACAS and Parameter file for ChIP-seq data Analysis 

## Test Data Sets

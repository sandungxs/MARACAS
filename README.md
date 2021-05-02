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


To install MARACAS follow the following steps:

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

3. Add to your .bashrc a new variable MARACAS indicating the path to the MARACAS scripts:  

```
echo "export MARACAS=$HOME/opt/MARACAS/scripts" >> .bashrc
source .bashrc
```

## How to use it?

To run our pipeline MARACAS you first have to create a parameter file. This parameter file is different depending if you want to analyze RNA-seq or ChIP-seq data. You can find examples of these files in the *test* folder. In the *test* folder you can also find a template for each type of parameter file. 

### Run MARACAS and Parameter file for RNA-seq data Analysis

To analyze RNA-seq data first access the folder 

### Run MARACAS and Parameter file for ChIP-seq data Analysis 

## Test Data Sets

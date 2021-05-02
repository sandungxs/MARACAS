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
* [Samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format](http://www.htslib.org/)
* [deepTools: tools for exploring deep sequencing data](https://deeptools.readthedocs.io/en/develop/index.html)
* [StringTie: Transcript assembly and quantification for RNA-Seq](https://ccb.jhu.edu/software/stringtie/)
* [R: free software environment for statistical computing and graphics with the following packages:](https://www.r-project.org/)
    * laslf
    * lasdjf



## How to use it?

## Test Data Sets

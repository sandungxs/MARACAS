
## Preprocess gff3 to generate gtf with gene_id and transcript_id 
vocar.gff3 <- read.table(file="volvox_carteri.gff3",header=F,quote = "#",as.is=T)
vocar.gtf <- vocar.gff3
head(vocar.gff3)

unique(vocar.gff3$V3)

for(i in 1:nrow(vocar.gff3))
{
  current.attributes <- strsplit(vocar.gff3$V9[i],split=";")[[1]]
  if(vocar.gff3$V3[i] == "gene")
  {
    gene.id <- strsplit(current.attributes[2],split="=")[[1]][2]
    vocar.gtf$V9[i] <- paste("gene_id", paste("\"",gene.id,"\";",sep=""))
  } else if(vocar.gff3$V3[i] == "mRNA")
  {
    gene.id <- substr(strsplit(current.attributes[5],split="=")[[1]][2],start = 1,stop = 15)
    transcript.id <- strsplit(current.attributes[2],"=")[[1]][2]
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  } else if(vocar.gff3$V3[i] == "exon")
  {
    gene.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 15)
    trancript.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 20)
    exon.number <- strsplit(strsplit(current.attributes[1],split="=")[[1]][2],split="exon.")[[1]][2]
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  }
  else if(vocar.gff3$V3[i] == "CDS")
  {
    gene.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 15)
    trancript.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 20)
    exon.number <- strsplit(strsplit(current.attributes[1],split="=")[[1]][2],split="CDS.")[[1]][2]
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  }
   else if(vocar.gff3$V3[i] == "five_prime_UTR" || vocar.gff3$V3[i] == "three_prime_UTR" || vocar.gff3$V3[i] == "CDS")
  {
    gene.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 15)
    transcript.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 20)
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  }
}

head(vocar.gtf)

write.table(x = vocar.gtf,file = "../../data/volvox_carteri/annotation/volvox_carteri.gtf",sep = "\t",row.names = F,col.names = F,quote = F)

## Generate TxDb package from gff3 file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

vocar.genome.data <- read.fasta(file = "../../data/dunaliella_salina/genome/dunaliella_salina.fa",seqtype = "DNA")
chromosome.names <- getName(vocar.genome.data)
chromosome.lengths <- sapply(X=getSequence(vocar.genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),value=c("https://phytozome.jgi.doe.gov/","v5.0"))

cre.txdb <- makeTxDbFromGFF(file = "../annotation/chlamydomonas_reinhardtii.gff3",format = "gff3",dataSource = "Phytozome",organism = "Chlamydomonas reinhardtii",taxonomyId = 3055,chrominfo = chromosomes.info,metadata = meta.data.info)

cre.txdb
genes(cre.txdb)

?makeTxDbPackage

makeTxDbPackage(txdb = cre.txdb, version = "0.1", maintainer = "Francisco J. Romero-Campero <fran@us.es>", author = "Francisco J. Romero-Campero")

install.packages("./TxDb.Creinhardtii.Phytozome/", repos=NULL)
## loading packages
library(ChIPseeker)
library(TxDb.Creinhardtii.Phytozome)
txdb <- TxDb.Creinhardtii.Phytozome
library(clusterProfiler)

files <- c("peaks_H3K27me3_before_starvation_1_peaks.narrowPeak","peaks_H3K27me3_sulfur_starvation_1_peaks.narrowPeak")
peak <- readPeakFile(files[[2]])
peak

covplot(peak, weightCol="X109")

covplot(peak, weightCol="X109",chrs="chromosome_1")
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-1000, 1000), color="red")

plotAvgProf(tagMatrix, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-1000, 1000), conf = 0.95, resample = 1000)

peakAnno <- annotatePeak(files[[2]], tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Creinhardtii.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

peak.annotation <- as.data.frame(peakAnno)
head(peak.annotation)
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-1000, 1000))

plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")

tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)
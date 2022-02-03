#working.directory <- arg[1]
working.directory <- "/home/ana/Documentos/GitHub/MARACAS/"
# control.condition <- args[2]
# experimental.condition <- args[3]
# fc.threshold <- as.numeric(args[4])
# q.val.threshold <- as.numeric(args[5])
# microalgae <- args[6]

#nuevo
#mapper <- args[7]
mapper <- "kallisto"
#num_samples <- args[8]
num_samples <- 2
setwd(working.directory)

gene.expression <- data.frame()
i<-1
for (i in 1:num_samples)
{
  current_sample <- paste(working.directory,"opt/kallisto_out_sample_", i, "/abundance.tsv", sep = "")
  #cambiar la ruta de opt a la carpeta de la muestra "samples/sample_", i, "/kallisto_out/abundance.tsv"
  kallisto_current_sample <- read.table(file="../opt/kallisto_out_sample1/abundance.tsv",header=T)
}
kallisto_sample_1 <- read.table(file="../opt/kallisto_out_sample1/abundance.tsv",header=T)
kallisto_sample_2 <- read.table(file="../opt/kallisto_out_sample2/abundance.tsv",header=T)

gene.expression <- data.frame(kallisto_sample_1[,1], kallisto_sample_1[,5], kallisto_sample_2[,5])
colnames(gene.expression) <- c("Transcripts ID", "sample_1", "sample_2")

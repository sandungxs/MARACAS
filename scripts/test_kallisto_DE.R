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

kallisto_sample1 <- read.table(file="../opt/kallisto_out_sample1/abundance.tsv",header=T)
kallisto_sample2 <- read.table(file="../opt/kallisto_out_sample2/abundance.tsv",header=T)

gene.expression <- data.frame(kallisto_sample1[,1], kallisto_sample1[,5], kallisto_sample2[,5])
colnames(gene.expression) <- c("Transcripts ID", "sample_1", "sample_2")

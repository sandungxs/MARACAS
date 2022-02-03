#working.directory <- arg[1]
working.directory <- "/home/ana/Documentos/GitHub/MARACAS/"
# control.condition <- args[2]
control.condition <- "no_iron"
# experimental.condition <- args[3]
experimental.condition <- "iron"
# fc.threshold <- as.numeric(args[4])
# q.val.threshold <- as.numeric(args[5])
# microalgae <- args[6]

#nuevo
#mapper <- args[7]
mapper <- "kallisto"

setwd(working.directory)

random_sample <- read.table(file=paste(working.directory,"opt/kallisto_out_sample_1/abundance.tsv",
                                       sep = ""), header=T)
gene.expression <- data.frame(random_sample[,1])
names <- c("transcript_ids")
for (i in 1:num_samples)
{
  current_sample <- paste(working.directory,"opt/kallisto_out_sample_", i, "/abundance.tsv", sep = "")
  #cambiar la ruta de opt a la carpeta de la muestra "samples/sample_", i, "/kallisto_out/abundance.tsv"
  kallisto_current_sample <- read.table(file=current_sample,header=T)
  gene.expression <- cbind(gene.expression,kallisto_current_sample[,5])
  names <- c(names, paste("sample_", i, sep=""))
}
colnames(gene.expression) <- names


## Authors: Francisco J. Romero-Campero
##          Ana Belén Romero-Losada
## Contact: Francisco J. Romero-Campero - email: fran@us.es

# working.directory <- "/home/fran/tmp/ostta_test/samples/"
# control.condition <- "iron"
# experimental.condition <- "no_iron"
# fc.threshold <- 2
# q.val.threshold <- 0.05

# working.directory <- "/home/fran/tmp/ngaditana/samples/"
# control.condition <- "high_N"
# experimental.condition <- "low_N"
# fc.threshold <- 2
# q.val.threshold <- 1


args <- commandArgs(trailingOnly=TRUE)

working.directory <- args[1]
control.condition <- args[2]
experimental.condition <- args[3]
fc.threshold <- as.numeric(args[4])
q.val.threshold <- as.numeric(args[5])

setwd(working.directory)

## Load libraries
library(ballgown)
library(genefilter)
library(FactoMineR)
library("factoextra")


# Load experimental design
experimental.design <- read.csv("experimental_design.csv",as.is=T)
experimental.design

## Sort samples 
sorted.samples <- sort(experimental.design$sample,ind=T)
indeces.sorted.samples <- sorted.samples$ix
experimental.design <- experimental.design[indeces.sorted.samples,]

number.samples <- nrow(experimental.design)
number.replicates <- table(experimental.design$condition)
control.replicates <- number.replicates[[control.condition]]
experimental.replicates <- number.replicates[[experimental.condition]]
sample.labels <- vector(mode="character",length = nrow(experimental.design))

control.indeces <- vector(mode="numeric",length=control.replicates)
experimental.indeces <- vector(mode="numeric",length=control.replicates)

j <- 1
k <- 1
for(i in 1:nrow(experimental.design))
{
  if(experimental.design$condition[i] == control.condition)
  {
    sample.labels[i] <- paste(control.condition,j,sep="_")
    control.indeces[j] <- i
    j <- j + 1
  } else if (experimental.design$condition[i] == experimental.condition)
  {
    sample.labels[i] <- paste(experimental.condition,k,sep="_")
    experimental.indeces[k] <- i
    k <- k + 1
  }
}

## Extract basic statistics for the result of each sample processing
## TODO!!!!!

## Load results from hisat2 + stringtie
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)

## Extract gene expression and name columns with sample.labels
gene.expression <- gexpr(bg.data)
colnames(gene.expression) <- sample.labels
head(gene.expression)
dim(gene.expression)

## Scatter plots 
number.samples <- nrow(experimental.design)
png(file="../results/scatter_plots.png",width = 1500,height = 1500)
par(mfrow=c(number.samples,number.samples))
for(i in 1:number.samples)
{
  for(j in 1:number.samples)
  {
    plot(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1),pch=19,cex=0.7,xlab=sample.labels[i],ylab=sample.labels[j],cex.lab=1.25)
    lines(x=c(-10,100),y=c(-10,100),col="red",lwd=2)
    text(x = 0, y =max(c(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1)))-2,pos = 4,
         paste0(round(x = 100*cor(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1)),digits=2),"%"),cex=2)
  }
}
dev.off()

## Boxplot before normalization
png(filename = "../results/boxplot_before_normalization.png")
boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,las=2,outline=F)
dev.off()

## PCA analysis
pca.gene.expression <- data.frame(colnames(gene.expression),t(gene.expression))
colnames(pca.gene.expression)[1] <- "condition"

res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
png(filename = "../results/pca_analysis_1.png")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70),main = "")
dev.off()

png(filename = "../results/pca_analysis_2.png")
fviz_pca_ind(res.pca, col.ind = experimental.design$condition, 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="Principal Components Analysis",
             show_legend=TRUE,show_guide=TRUE) 
dev.off()

## Hierarchical clustering
#res.hcpc <- HCPC(res.pca, graph=FALSE)    

#png(filename = "../results/hierarchical_clustering.png")
#plot(res.hcpc, choice ="tree")
#dev.off()

## Apply upper quantile normalization
upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

mean.upper.quantiles <- mean(upper.quantiles)

for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- (gene.expression[,i] / upper.quantiles[i]) * mean.upper.quantiles
}

## Log2 transformation
log.gene.expression <- log2(gene.expression+1)

png(filename = "../results/boxplot_after_normalization.png")
boxplot(log.gene.expression,col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,outline=F,las=2,main="Normalized Gene Expression")
dev.off()

## Alternativamente la normalización se puede realizar con el paquete normalyzer. 
# ## En este punto ballgown no realiza ninguna normalización de los datos. 
# ## Utilizamos el paquete de R Normalyzer para esta tarea. Para ello es neceario generar
# ## un fichero con un formato específico.
# 
# normalyzer.data <- data.frame(rownames(gene.expression),gene.expression)
# dim(normalyzer.data)
# head(normalyzer.data)
# 
# colnames(normalyzer.data) <- NULL
# rownames(normalyzer.data) <- NULL
# 
# normalyzer.table <- rbind(c(0,rep(1:2,each=2)),                                  
#                     rbind(c("Gene",colnames(gene.expression)),normalyzer.data))
# 
# head(normalyzer.table)
# 
# write.table(normalyzer.table,file="normalyzer_table.tab",col.names = F, quote = F, row.names = F, sep="\t")
# 
# library(Normalyzer)
# library(grid)
# normalyzer(datafile = "normalyzer_table.tab", getjob = "data_normalization")
# 
# normalized.data <- read.table(file="data_normalization/Quantile-normalized.txt", header=T)
# head(normalized.data)
# 
# normalized.data[is.na(normalized.data)] <- 0
# 
# ## Testeamos la normalización
# boxplot(normalized.data[,2:5],col=rep(c("red","blue"),each=2))

## Compute the mean expression matrix
if(length(control.indeces) > 1)
{
  control <- rowMeans(log.gene.expression[,control.indeces])  
} else
{
  control <- as.vector(log.gene.expression[,control.indeces])
  names(control) <- rownames(log.gene.expression)
}

if(length(experimental.indeces) > 1)
{
  experimental <- rowMeans(log.gene.expression[,experimental.indeces])
} else
{
  experimental <- as.vector(log.gene.expression[,experimental.indeces])
  names(experimental) <- rownames(log.gene.expression)
}


mean.expression <- matrix(c(control,experimental),ncol=2)
colnames(mean.expression) <- c(control.condition,experimental.condition)
rownames(mean.expression) <- names(control)
head(mean.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
plot(control,experimental,pch=19,cex=0.7,xlab=control.condition,ylab=experimental.condition,cex.lab=1.25)

##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

## Specification of the experimental design
factor.experimental.design <- vector(mode="numeric",length=nrow(experimental.design))
factor.experimental.design[control.indeces] <- 1
factor.experimental.design[experimental.indeces] <- 2

limma.experimental.design <- model.matrix(~ -1+factor(factor.experimental.design))
colnames(limma.experimental.design) <- c(control.condition,experimental.condition)

## Linear model fit
linear.fit <- lmFit(log.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(contrasts = paste(experimental.condition,control.condition,sep="-"),
                                 levels=c(control.condition,experimental.condition))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract results
de.results <- topTable(contrast.results, number=7507,coef=1,sort.by="logFC")
head(de.results)

fold.change <- de.results$logFC
q.values <- de.results$adj.P.Val
genes.ids <- rownames(de.results)

names(fold.change) <- genes.ids
names(q.values) <- genes.ids

activated.genes <- genes.ids[fold.change > log2(fc.threshold) & q.values < q.val.threshold]
repressed.genes <- genes.ids[fold.change < - log2(fc.threshold) & q.values < q.val.threshold]

length(activated.genes)
length(repressed.genes)

write(x = activated.genes, file = "../results/activated_genes.txt")
write(x = repressed.genes, file = "../results/repressed_genes.txt")

png(filename = "../results/scatter_plot_control_vs_experimental.png")
plot(control,experimental,pch=19,cex=0.7,col="grey",xlab=control.condition,ylab=experimental.condition,cex.lab=1.25)
points(control[activated.genes],experimental[activated.genes],pch=19,cex=0.7,col="red")
points(control[repressed.genes],experimental[repressed.genes],pch=19,cex=0.7,col="blue")
dev.off()

log10.qval <- -log10(q.values)

png(filename = "../results/volcano_plot.png")
plot(fold.change,log10.qval,pch=19,cex=0.7,col="grey", xlab="Fold Change", ylab="-log10(q-value)",cex.lab=1.5)
points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)
dev.off()

## Código para desarrollar una función gráfico de barras
gene <- activated.genes[2]

original.data <- 2^log.gene.expression - 1

control.expr.vals <- unlist(c(original.data[gene, control.indeces]))
experimental.expr.vals <- unlist(c(original.data[gene, experimental.indeces]))

mean.control <- mean(control.expr.vals)
mean.experimental <- mean(experimental.expr.vals)

if(length(control.indeces) > 1)
{
  sd.control <- sd(control.expr.vals)  
} else
{
  sd.control <- 0
}
  
if(length(experimental.indeces) > 1)
{
  sd.experimental <- sd(experimental.expr.vals)  
} else
{
  sd.experimental <- 0
}

means <- c(mean.control, mean.experimental)
sds <- c(sd.control, sd.experimental)

arrow.top <- means + sds
arrow.bottom <- means - sds

png(filename = "../results/barplot.png")
xpos <- barplot(means,ylim=c(0,1.5*max(arrow.top)),col=c("red","blue"),names.arg = c(control.condition,experimental.condition),ylab="FPKM",cex.lab=1.5,main=gene,cex.main=2)
arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.1,lwd=1.5)
points(rep(xpos[1],length(control.expr.vals))+0.1,control.expr.vals)
points(rep(xpos[2],length(experimental.expr.vals))+0.1,experimental.expr.vals)
dev.off()

       

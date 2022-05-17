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

# working.directory <- "/home/fran/tmp/local_test/samples/"
# control.condition <- "iron"
# experimental.condition <- "no_iron"
# fc.threshold <- 2
# q.val.threshold <- 1
# microalgae <- "Micromonas pusilla CCMP1545"

# working.directory <- "/home/fran/tmp/rna_seq_test/samples/"
# control.condition <- "iron"
# experimental.condition <- "no_iron"
# fc.threshold <- 2
# q.val.threshold <- 1
# microalgae <- "test"
# mapper <- "kallisto"


args <- commandArgs(trailingOnly=TRUE)

working.directory <- args[1]
control.condition <- args[2]
experimental.condition <- args[3]
fc.threshold <- as.numeric(args[4])
q.val.threshold <- as.numeric(args[5])
microalgae <- args[6]
mapper <- args[7]
setwd(working.directory)

## Load libraries
library(ballgown)
library(genefilter)
library(FactoMineR)
library("factoextra")
library(ggplot2)
# library(MetBrewer)
library(cowplot)
# library(plotly)

# Load experimental design
experimental.design <- read.csv("experimental_design.csv",as.is=T)
experimental.design
number.samples <- nrow(experimental.design)

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

if (mapper == "hisat2")
{
  ## Load results from hisat2 + stringtie
  bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)
  
  ## Extract gene expression and name columns with sample.labels
  gene.expression <- gexpr(bg.data)
  colnames(gene.expression) <- sample.labels
  head(gene.expression)
  dim(gene.expression)
  df.gene.expression <- data.frame(row.names(gene.expression),gene.expression)
  colnames(df.gene.expression)[1] <- "geneID"
  rownames(df.gene.expression) <- NULL
  
  write.table(x = df.gene.expression,file = "../results/gene_expression.tsv",
              quote = F,sep = "\t",row.names = F)
} else if (mapper == "kallisto")
{
  #Extract transcripts abundance from each kallisto output to create a gene expression table
  random_sample <- read.table(file="./sample_1/kallisto_out/abundance.tsv", header=T)
  gene.expression <- data.frame(matrix(NA, nrow=length(random_sample[,1]), ncol=number.samples))
  tpm.expression <- data.frame(matrix(NA, nrow=length(random_sample[,1]), ncol=number.samples))
  gene.count <- data.frame(matrix(NA, nrow=length(random_sample[,1]), ncol=number.samples))
  names <- c()
  
  for (i in 1:number.samples)
  {
    kallisto_current_sample <- read.table(file=paste(c("./sample_", i,"/kallisto_out/abundance.tsv"),collapse=""),header=T)
    gene.expression[,i] <- kallisto_current_sample[,5]
    tpm.expression[,i] <- kallisto_current_sample[,5]
    gene.count[,i] <- kallisto_current_sample[,4]
    names <- c(names, paste("sample_", i, sep=""))
  }
  colnames(gene.expression) <- names
  rownames(gene.expression) <- random_sample[,1]
  write.table(gene.expression, file = "../results/gene_expression.tsv", quote = F )
  
  colnames(tpm.expression) <- names
  rownames(tpm.expression) <- random_sample[,1]
  write.table(tpm.expression, file = "../results/transcript_count_matrix.csv", quote = F )
  
  colnames(gene.count) <- names
  rownames(gene.count) <- random_sample[,1]
  write.table(gene.count, file = "../results/gene_count_matrix.csv", quote = F )
}


## Scatter plots 

scatterplot_function<-function(x,y)
{
  print(ggplot(as.data.frame(log2(gene.expression)+1), 
         aes(x=log2(gene.expression[,x])+1, y=log2(gene.expression[,y])+1)) + 
    #geom_point(size=1, color=met.brewer("Cassatt2",8) [6],
    geom_point(size=1, color="#7fa074")+
               #aes(text= paste0("</br> X: ",round(log2(gene.expression[,x])+1,digits = 5),
               #                 "</br> Y: ", round(log2(gene.expression[,y])+1,digits = 5)))) +
    xlim(0,NA) +
    ylim(0,NA) +
    geom_smooth(method=lm,color="darkgreen",size=0.7) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          axis.line.x.bottom = element_line(color = 'black'),
          axis.line.y.left   = element_line(color = 'black'),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 13)) +
    xlab(sample.labels[x]) +
    ylab(sample.labels[y]) +
    annotate(geom = "text", x = 1.5, y = max(c(log2(gene.expression[,x]+1),log2(gene.expression[,y]+1)))-2, 
             label = paste(c(round(100*cor(gene.expression[,x],
                                           gene.expression[,y]),
                                   digits = 2),
                             "%"), collapse=""),size=4))
}

png(file="../results/scatter_plots.png",width = 1500,height = 1500)
#par(mfrow=c(number.samples,number.samples))

#for(i in 1:number.samples)
#{
#  for(j in 1:number.samples)
#  {
#    plot(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1),pch=19,
#         cex=0.7,xlab=sample.labels[i],ylab=sample.labels[j],cex.lab=1.5)
#    lines(x=c(-10,100),y=c(-10,100),col="red",lwd=2)
#    text(x = 0, y =max(c(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1)))-2,pos = 4,
#         paste0(round(x = 100*cor(log2(gene.expression[,i]+1),log2(gene.expression[,j]+1)),digits=2),"%"),cex=2)
#  }
#}

scatterplot_function(2,4)

dev.off()

# myPlots<- list()

# for(i in 1:number.samples)
#  {
#    for(j in 1:number.samples)
#    {
#      myPlots<-c(myPlots,list())
#    }
#}

#print(ggarrange(plotlist = myPlots, nrow = number.samples, ncol = number.samples))

#dev.off()

## Boxplot before normalization
if(mapper == "hisat2")
{
  png(filename = "../results/boxplot_before_normalization.png")
# boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,las=2,outline=F)
  print(ggplot(stack(as.data.frame(log2(gene.expression+1))), aes(x = ind, y = values, fill = ind)) +
        stat_boxplot(geom = 'errorbar') +
        geom_boxplot(outlier.shape = NA) + 
        scale_y_continuous(limits = quantile(log2(gene.expression[,1]+1), c(0.1, 0.9))) +
        # scale_fill_manual(values=c(met.brewer("Klimt",number.samples))) +
        scale_fill_manual(values=rainbow(ncol(gene.expression))) +
        theme(legend.position="none", axis.text.x = element_text(angle = 90),
              # plot.title = element_text(hjust = 0.5, size = 13),
              panel.background = element_rect(fill = "white"),
              axis.title = element_text(size = 11),
              axis.text = element_text(size=9),
              panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"),
              axis.line.x.bottom = element_line(color = 'black'),
              axis.line.y.left   = element_line(color = 'black'),
              panel.border = element_blank()
        )  +
        xlab("Samples") +
        ylab("Gene Expression"))

dev.off()
}else if (mapper == "kallisto")
{
  png(filename = "../results/boxplot_before_normalization.png")
  # boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(TPM + 1)",cex.lab=1.5,las=2,outline=F)
  print(ggplot(stack(as.data.frame(log2(gene.expression+1))), aes(x = ind, y = values, fill = ind)) +
          stat_boxplot(geom = 'errorbar') +
          geom_boxplot(outlier.shape = NA) + 
          scale_y_continuous(limits = quantile(log2(gene.expression[,1]+1), c(0.1, 0.9))) +
          scale_fill_manual(values=rainbow(ncol(gene.expression))) +
          # scale_fill_manual(values=c(met.brewer("Klimt",number.samples))) +
          theme(legend.position="none", axis.text.x = element_text(angle = 90), 
                plot.title = element_text(hjust = 0.5, size = 13), 
                panel.background = element_rect(fill = "white"),
                axis.title = element_text(size = 11),
                axis.text = element_text(size=9),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                axis.line.x.bottom = element_line(color = 'black'),
                axis.line.y.left   = element_line(color = 'black'),
                panel.border = element_blank()
          )  +
          xlab("Samples") +
          ylab("TPM"))
  dev.off()
}


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

if (mapper == "hisat2")
{
   png(filename = "../results/boxplot_after_normalization.png")
   # boxplot(log.gene.expression,col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,outline=F,las=2,main="Normalized Gene Expression")
   print(ggplot(stack(as.data.frame(log.gene.expression)), aes(x = ind, y = values, fill = ind)) +
     stat_boxplot(geom = 'errorbar') +
     geom_boxplot(outlier.shape = NA) + 
     scale_y_continuous(limits = quantile(log.gene.expression[,1], c(0.1, 0.9))) +
       scale_fill_manual(values=rainbow(ncol(gene.expression))) +
     # scale_fill_manual(values=c(met.brewer("Klimt",number.samples))) +
     theme(legend.position="none", axis.text.x = element_text(angle = 90), 
           plot.title = element_text(hjust = 0.5, size = 13), 
           panel.background = element_rect(fill = "white"),
           axis.title = element_text(size = 11),
           axis.text = element_text(size=9),
           panel.grid.major = element_line(colour = "white"),
           panel.grid.minor = element_line(colour = "white"),
           axis.line.x.bottom = element_line(color = 'black'),
           axis.line.y.left   = element_line(color = 'black'),
           panel.border = element_blank()
     ) +
     xlab("Samples") +
     ylab("log2(FPKM + 1)"))
   dev.off()
} else if (mapper == "kallisto")
{
   png(filename = "../results/boxplot_after_normalization.png")
   # boxplot(log.gene.expression,col=rainbow(ncol(gene.expression)),ylab="log2(TPM + 1)",cex.lab=1.5,outline=F,las=2,main="Normalized Gene Expression")
   print(ggplot(stack(as.data.frame(log.gene.expression)), aes(x = ind, y = values, fill = ind)) +
           stat_boxplot(geom = 'errorbar') +
           geom_boxplot(outlier.shape = NA) + 
           scale_y_continuous(limits = quantile(log.gene.expression[,1], c(0.1, 0.9))) +
           scale_fill_manual(values=rainbow(ncol(gene.expression))) +
           # scale_fill_manual(values=c(met.brewer("Klimt",number.samples))) +
           theme(legend.position="none", axis.text.x = element_text(angle = 90), 
                 plot.title = element_text(hjust = 0.5, size = 13), 
                 panel.background = element_rect(fill = "white"),
                 axis.title = element_text(size = 11),
                 axis.text = element_text(size=9),
                 panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 axis.line.x.bottom = element_line(color = 'black'),
                 axis.line.y.left   = element_line(color = 'black'),
                 panel.border = element_blank()
           ) +
           xlab("Samples") +
           ylab("log2(TPM + 1)"))
   dev.off()
}

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
#head(mean.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
# plot(control,experimental,pch=19,cex=0.7,xlab=control.condition,ylab=experimental.condition,cex.lab=1.25)

# print(ggplot(mean.expression,aes(x=control.condition,y=experimental.condition)) + 
#  geom_point(size=1,colour= "#924099") +
  # geom_point(size=1,colour=c(met.brewer("Klimt",1))) +
#  xlab(control.condition) +
#  ylab(experimental.condition) +
#  theme(axis.title = element_text(size = 15)))


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
#head(de.results)

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

# Scatterplot DEGs
png(filename = "../results/scatter_plot_control_vs_experimental.png")
#  plot(control,experimental,pch=19,cex=0.7,col="grey",xlab=control.condition,ylab=experimental.condition,cex.lab=1.25)
#  points(control[activated.genes],experimental[activated.genes],pch=19,cex=0.7,col="red")
#  points(control[repressed.genes],experimental[repressed.genes],pch=19,cex=0.7,col="blue")

mean.expression<-as.data.frame(mean.expression)
# head(mean.expression)
mean.expression[,"gene_type"] <- "ns" 

mean.expression[which(mean.expression$high_light - mean.expression$control > fc.threshold),"gene_type"] <- "activado"
mean.expression[which(mean.expression$high_light - mean.expression$control < -(fc.threshold)),"gene_type"] <- "reprimido"

volcol <- c("#dd5129", "#0f7ba2","grey33")
names(volcol) <- c("activado","reprimido","ns")

print(ggplot(mean.expression, aes(x=control, y=experimental, color=gene_type)) + 
  geom_point() +
  scale_color_manual(values=volcol) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5, size = 13), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size=9),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank()
  ))

dev.off()

log10.qval <- -log10(q.values)

png(filename = "../results/volcano_plot.png")
# plot(fold.change,log10.qval,pch=19,cex=0.7,col="grey", xlab="Fold Change", ylab="-log10(q-value)",cex.lab=1.5)
# points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
# points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)

de.results[,"gene_type"] <- "ns"
de.results[,"gene_name"] <- rownames(de.results)
de.results[which(de.results$logFC > fc.threshold & de.results$adj.P.Val < q.val.threshold),"gene_type"] <- "activado"
de.results[which(de.results$logFC< -fc.threshold & de.results$adj.P.Val < q.val.threshold),"gene_type"] <- "reprimido"

 volcol <- c("#dd5129", "#0f7ba2","grey33")
# volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
 names(volcol) <- c("activado","reprimido","ns")

print(ggplot(as.data.frame(de.results), aes(x=logFC, y=-log10(adj.P.Val),
                                                   color=gene_type)) +
  geom_point(size = 1) +
   scale_colour_manual(values = volcol) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19)))

dev.off()

## Código para desarrollar una función gráfico de barras
if(length(activated.genes) > 2)
{
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
}

## Generation of Rmd file for final report
output.file <- "../results/DE_report.Rmd"

## Document header
header <- paste0(paste0("# **Differential Gene Expression Report for ",
                        microalgae),
                 "**")

## Introduction test
write(x = header,file = output.file,append = F)
write(x = "## Experimental Design and Global Statistics",file = output.file,append = T)

intro.line <- paste(c("This report was automatically generated by **MARACAS 
                      (MicroAlgae RNA-seq and ChIP-seq Analysis)** as a 
                      differential expression analysis for the microalgae **",
                      microalgae, "**. The experimental design considers as control
                      condition **", control.condition, "** with ",
                      control.replicates, " replicates and as experimental condition **",
                      experimental.condition, "** with ", experimental.replicates,
                      " replicates. A summary of the processing of each individual
                      sample is presented below. Click on the corresponding links to 
                      download a **quality control analysis** or a file in **bigWig**
                      format containing information about the **mapping signal** for 
                      each sample:\n"), collapse="")
write(x = intro.line,file = output.file,append = T)

## Experimental design
number.samples <- nrow(experimental.design)

for(i in 1:number.samples)
{
  write(x = paste(c("* **", experimental.design$sample[i], "** for condition ",
                    "**", experimental.design$condition[i], "**:"),collapse=""),file = output.file,
        append = T)
  write(x=paste(c("[Quality Control analysis](../samples/",
                  experimental.design$sample[i],
                  "/sample_1_fastqc.html), "), collapse=""), 
        file=output.file, append=T)
  if( mapper == "hisat2")
  {
    write(x=paste(c("[BigWig file with mapping signal](../samples/",
                    experimental.design$sample[i],"/sample.bw)"),
                  collapse=""), 
          file=output.file, append=T)
  }
  
  if( mapper == "hisat2" )
  {
    
    mapping.stats <- readLines(con = paste(c( experimental.design$sample[i],
                                              "/mapping_stats"),collapse=""),n = 6) 
    write(x = "\n",file=output.file, append=T)
    for(j in 1:6)
    {
      write(x=paste(c("\t",mapping.stats[j],"\n"),collapse = ""),file=output.file, append=T)
    }
  }
  #write(x = paste(c("\t",mapping.stats[6]),collapse = ""),file=output.file, append=T)
  write(x = "\n",file=output.file, append=T)
}

write(x = "\n",file=output.file, append=T)
write(x = "\n",file=output.file, append=T)

if (mapper == "hisat2")
{
   write(x="[**Click here to download a matrix in tab-separated value format containing 
   estimates for gene expression computed from your RNA-seq data measured as FPKM. Rows represent genes
      and columns conditions.**](./gene_expression.tsv)", file=output.file, append=T)
   write(x = "\n",file=output.file, append=T)
}

write(x="[**Click here to download a matrix in tab-separated value format containing 
estimates for gene expression computed from your RNA-seq data measured as TPM. Rows represent genes
      and columns conditions.**](./transcript_count_matrix.csv)", file=output.file, append=T)
write(x = "\n",file=output.file, append=T)
write(x="[**Click here to download a matrix in tab-separated value format containing 
estimates for gene expression computed from your RNA-seq data measured as raw read count. Rows represent genes
      and columns conditions.**](./gene_count_matrix.csv)", file=output.file, append=T)
write(x = "\n",file=output.file, append=T)

## Global Gene Expression
write(x = "\n", file=output.file, append=T)
write(x="The barplot below represents the upper quantile normalized 
      global gene expression distibrution for each sample. The global gene
      expression distribution should be comparable among all samples:\n", file=output.file, append=T)

write(x = "<center>", file=output.file, append=T)
write(x = "![Normalized Global Gene Expression](./boxplot_after_normalization.png){ width=50% }", file=output.file, append=T)
write(x = "</center>", file=output.file, append=T)

## Principal components analysis
write(x = "\n", file=output.file, append=T)
write(x="A graphical representation of a **Principal Component Analysis (PCA)** 
      is presented below. Control and experimental samples should cluster into 
      different and distinct groups:\n", file=output.file, append=T)

write(x = "<center>", file=output.file, append=T)
write(x = "![Principal Components Analysis](./pca_analysis_2.png){ width=50% }", file=output.file, append=T)
write(x = "</center>", file=output.file, append=T)

## Samples Scatter plots.
write(x="Scatter plots comparing gene expression between samples is presented
      below. Samples constituting replicates for the same condition should have
      a strong diagonal signal and high correlation value.\n", file=output.file, append=T)

write(x = "<center>", file=output.file, append=T)
write(x = "![Sample Scatter Plots](./scatter_plots.png){ width=100% }", file=output.file, append=T)
write(x = "</center>", file=output.file, append=T)

write(x = "## Differential Gene Expression Analysis",file = output.file,append = T)

write(x = paste(c("In this section we present the results from a differential gene expression 
                   analysis performed using a **q-value of ", q.val.threshold, "** and 
                   a **fold-change of ", fc.threshold, "**. The scatter and volcano
                   plots below represent the gene differential expression. In both
                  graphs the red and blue dots represent activated and repressed
                  genes respectively. The grey dots represent genes that not 
                  differentially expressed. In the **scatter plot**, the x-coordinate
                  represents the log2 gene expression in the control condition and 
                  the y-coordinate represents the log2 gene expression in the 
                  experimental condition. Grey dots in the diagonal stand for non
                  differentially expressed genes, red dots in the upper triangle
                  represent activated genes and blue dots inthe lower triangle
                  specify repressed genes. In the **volcano plot**, the x-coordinate
                  represents the **log2 gene expression fold-change** and the 
                  y-coordinate represent the **minus log10 q-value**. Similarly,
                  grey, red and blue dots represent non differentially expressed,
                  activated and repressed genes respectively.\n 
                  " ),collapse=""),
      file = output.file,append = T)

write(x = "<center>", file=output.file, append=T)
write(x = "![Scatter Plot](./scatter_plot_control_vs_experimental.png){ width=50% }", file=output.file, append=T)
write(x = "</center>", file=output.file, append=T)

write(x = "<center>", file=output.file, append=T)
write(x = "![Scatter Plot](./volcano_plot.png){ width=50% }", file=output.file, append=T)
write(x = "</center>", file=output.file, append=T)

number.activated.genes <- length(activated.genes)
number.repressed.genes <- length(repressed.genes)

write(x = paste(c("Specifically, we detected **", number.activated.genes, " activated
                  genes** and **", number.repressed.genes, " repressed genes**. Click
                  on the links below to download text files containing the gene IDs for
                  activated and repressed genes. These files can be used to perform
                  a functional enrichment analysis using our tool **ALGAEFUN (
                  Microalge Functional Annotation Enrichment Analysis)**: \n
                  "),collapse=""),
      file = output.file,append = T)

write(x = "* [**Click here to download the list of activated genes.**](./activated_genes.txt)\n",
      file = output.file,append = T)
write(x = "* [**Click here to download the list of repressed genes.**](./repressed_genes.txt)\n",
      file = output.file,append = T)

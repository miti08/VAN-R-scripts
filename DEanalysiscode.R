##setwd to the folder containing the dataset
##composer:HoangAnhLeP
##Institute: Brain Health Lab - International University, Vietnam
##composingdate:19-03-2023

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(edgeR)
library(DESeq2)
library(ggplot2)
library(limma)
library(gplots)
library(AnnotationDbi)

library(Glimma)
library(RColorBrewer)
##loaddata and create DGE list
data <- read.table("merged_sumrawcounts.txt",sep=",",header=TRUE) ##input the rawcounts data post-cleaning
dge <- DGEList(data,group=substr(colnames(data),0,1))

##Extract sample info
Counts <- read.table("merged_sumrawcounts.txt",sep=",",header=TRUE)
Counts <- Counts[which(rowSums(Counts) > 0),] #omit the genes without counts
Counts
condition <- factor(c("C","C","C","C","C","C","C","C","C","C","T","T","T","T","T","T","T","T","T","T")) #label subjects' condition
coldata <- data.frame(row.names = colnames(Counts), condition)
coldata

##run DE Analysis
dds <- DESeqDataSetFromMatrix(countData = dge$counts, colData = coldata, design = ~condition)
dds1 <- estimateSizeFactors(dds) #estimate the size factors of each subject
dds <- DESeq(dds)
sizeFactors(dds1) #view the size factors dedicated to each subject
normalized_counts <- counts(dds1, normalized=TRUE) #run DESeq2 normalization on the rawcounts
normalized_counts <- normalized_counts[order(rowSums(normalized_counts)),] #sort the data.frame by total reads per gene
write.csv(normalized_counts, file="Counts_significantresults_reorder.csv")

##Get variance of dataset
vsdata <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
plotPCA(vsdata, intgroup = "condition")
plotDispEsts(dds,log="xy",legend=TRUE,cex=0.5)

##Report results from DE analysis
res <- results(dds, contrast = c("condition", "T", "C"))
res
sigs <- na.omit(res) ##omit na points
write.csv(sigs,file="DESeq_result_unnormalize.csv")

##select significant data
library(dplyr)
sigss <- sigs[sigs$padj<0.05,]
sigss
write.csv(sigss,file="DESeq_significantresults_unnormalized.csv")

##Clasify into up and down regulated genes
dess <- c()
Ten <- rownames(sigss)
for (i in (1:length(Ten)))
{if (sigss$log2FoldChange[i]>0)
  {dess <- c(dess,"UP")} else {dess <- c(dess,"DOWN")}}
df <- data.frame(Ten, sigss$log2FoldChange, sigss$padj, dess)
write.csv(df,file="DESeq_expressionlevel_unnormalized.csv")

ls
df <- as.data.frame(sigs)

##draw volcano plot from log2foldchange and p-value
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
VP <- data.frame(sigs@rownames,sigs$log2FoldChange,sigs$padj) #create data.frame contains log2Foldchance and p-val for volcanoplot
rownames(VP)<-NULL #remove numeric labels of rows
VP
VPlot<-EnhancedVolcano(toptable = VP,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "sigs.log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "sigs.padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(VP),labSize=5,ylim=c(0,50),xlim=c(-30,30))

##cut dataframe choosing significant gene counts only
library(dplyr)
countsel<- Counts[rownames(Counts) %in% Ten, ]


##Seperate the counts into up and down
K <-read.csv("DESeq_expressionlevel_unnormalized.csv")
library(dplyr)
uplist <- K[K$dess %in% "UP",]           #select the genes that upregulate
uplist1 <- uplist$Ten                    #build the list of regulated genes
downlist <- K[K$dess %in% "DOWN",]       #select the genes that downregulate
downlist1 <-downlist$Ten                 #build the list of downregulated genes
countsel <- read.csv("Counts_significantresults_reorder.csv") 
countup <- countsel[countsel$X %in% uplist1,]                  #select normalized counts data belong to the upregulated list
countdown <- countsel[countsel$X %in% downlist1,]              #select normalized counts data belong to the downregulated list
write.csv(countup,file="upregulated_genecounts_normalized.csv")
write.csv(countdown,file="downregulated_genecounts_normalized.csv")

##get percentile of dataframe to build color scale for heatmap
library(rstatix)
datacountup <- countup %>% select(colnames(countup),-contains("X")) #remove the genename column
datacountup <- as.vector.data.frame(datacountup) #convert the dataframe into vector
uppercentile <- quantile(unlist(datacountup), probs = c(0.2,0.4,0.5,0.6,0.7,0.8,0.9)) #unlist the vector to get atomic elements
uppercentile #print the results
datacountdown <- countdown %>% select(colnames(countdown),-contains("X"))
datacountdown <- as.vector.data.frame(datacountdown)
downpercentile <- quantile(unlist(datacountdown), probs = c(0.2,0.4,0.5,0.6,0.7,0.8,0.9))
downpercentile

res_lfc <- subset(sigss, abs(log2FoldChange) > 1) 
topVarGenes <- head(order(rowVars(assay(vsdata), useNames = T), decreasing = T), 50)
mat  <- assay(vsdata)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
top20DEGs <- order(abs(res_lfc$log2FoldChange), decreasing=TRUE)


coldata2<-data.frame(sample=colnames(Counts), condition=c("control","control","control","control","control","control","control","control","control","control","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer","alzheimer"))
annot_col <- coldata2 %>%
  column_to_rownames('sample') %>%
  dplyr::select(condition) %>%
  as.data.frame()

pheatmap(mat, cluster_rows=TRUE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=annot_col, fontsize_row=6)

### WGCNA script

## Required libraries
library(WGCNA)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(ggpubr)
#EXPRESSION MATRIX INPUT---------------------------------
data <- read.table("merged_sumrawcounts.txt",sep=",",header=TRUE)
condition <- factor(c("C","C","C","C","C","C","C","C","C","C","T","T","T","T","T","T","T","T","T","T")) #label subjects' condition
coldata <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~condition)
dds <- DESeq(dds)
vsdata <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
expr_data_vst<-t(assay(vsdata))

options(stringsAsFactors = FALSE);
phenotype_mat <- read.csv("metadata.csv", header=TRUE,row.names=1)
datExpr<-as.data.frame(expr_data_vst)
#Get list of upregulated and downregulated genes from DE analysis
DEgenes <- read.csv("DESeq_expressionlevel_unnormalized.csv", header=TRUE)
upgenes<-DEgenes %>% filter(dess=="UP") %>% select("genes")
downgenes<-DEgenes %>% filter(dess=="DOWN") %>% select("genes")

## NETWORK CONSTRUCTION----------------------------------------------

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  datExpr,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5,
  networkType = 'signed',
  # corFnc = bicor,
  # corOptions = list(maxPOutliers=0.05)
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 12 #choose a power level that is close to the threshold
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)        
netwk <- blockwiseModules(datExpr,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 20000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

### MODULE ANALYSES-----------------------------
load("D:/VAN_data/ER-block.1.RData")

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_eigengenes <- netwk$MEs

# Print out a preview
head(module_eigengenes)

#number of genes in each module - plot

yellowGenes<-names(datExpr)[mergedColors=="yellow"]
turqGenes<-names(datExpr)[mergedColors=="turquoise"]
brownGenes<-names(datExpr)[mergedColors=="brown"]

moduleGeneNumbers <- as.data.frame(table(mergedColors))
barplot(moduleGeneNumbers$Freq, 
        names.arg = moduleGeneNumbers$mergedColors,
        ylab="Number of genes",
        main="Number of genes in each module",
        ylim=c(0,150),
        col=c("blue","brown","green","grey","red","turquoise","yellow"))

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels 
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0) 
moduleTraitCor = cor(MEs, phenotype_mat, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("AD Status","Age","Sex"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


### INTRAMODULAR ANALYSES-----------------------------
# Define variable weight containing the weight column of datTrait
alz = as.data.frame(phenotype_mat$alzheimer);
names(alz) = "alzheimer"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, alz, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(alz), sep="");
names(GSPvalue) = paste("p.GS.", names(alz), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = mergedColors==module;

MG_x<-abs(geneModuleMembership[moduleGenes, column])
MG_y<-abs(geneTraitSignificance[moduleGenes, 1])
MMvsGS<-tibble(label=names(datExpr)[moduleGenes],x=MG_x,y=MG_y)
getMMGS<-function(module,...){
  column = match(module, modNames);
  moduleGenes = mergedColors==module;
  MG_x<-abs(geneModuleMembership[moduleGenes, column])
  MG_y<-abs(geneTraitSignificance[moduleGenes, 1])
  MMvsGS<-tibble(label=names(datExpr)[moduleGenes],x=MG_x,y=MG_y)
  return(MMvsGS)
}
MMGS_turq<-getMMGS("turquoise")
MMGS_brown<-getMMGS("brown")
MMGS_yellow<-getMMGS("yellow")
bigMMlist<-MMGS_turq %>% full_join(MMGS_brown) %>% full_join(MMGS_yellow)

MMgenes<-filter(MMvsGS, x>0.8)
DEGs<-append(downgenes$genes,upgenes$genes)

venn_turq<-ggVennDiagram(list(MMgenes_turq$label,DEGs),
              category.names=c("Genes with MM>0.8","DEGs"),
              show_intersect=TRUE,set_size=1,label_size=1)+labs(title="Turquoise module")+scale_fill_distiller(palette = "RdBu")
venn_brown<-ggVennDiagram(list(MMgenes_brown$label,DEGs),
                         category.names=c("Genes with MM>0.8","DEGs"),
                         show_intersect=F)+labs(title="Brown module")+scale_fill_distiller(palette = "RdBu")
venn_yellow<-ggVennDiagram(list(MMgenes_yellow$label,DEGs),
                         category.names=c("Genes with MM>0.8","DEGs"),
                         show_intersect=F)+labs(title="Yellow module")+scale_fill_distiller(palette = "RdBu")

MMgenes_turq<-filter(MMGS_turq, x>0.8)
MMgenes_brown<-filter(MMGS_brown, x>0.8)
MMgenes_yellow<-filter(MMGS_yellow, x>0.8)
HubGenes_yellow<-MMgenes_yellow %>% filter(label %in% upgenes$genes)


AllHubs<-c(HubGenes_brown$label,HubGenes_turq$label,HubGenes_yellow$label)
AllHubs_MM<-c(MMgenes_turq$label,MMgenes_brown$label,MMgenes_yellow$label)
##LABEL WITH UPREGULATED AND DOWNREGULATED GENES
downoverlap<-MMvsGS %>% filter(label %in% downgenes$genes)
upoverlap <- MMvsGS %>% filter(label %in% upgenes$genes)
p1<-ggplot(MMGS_turq, aes(x=x, y=y,label=label)) + 
  geom_point(shape="circle",colour="black",size=2,alpha=0.5) +
  geom_point(data=MMGS_turq %>% filter(label %in% downgenes$genes),aes(x=x,y=y),color="red") +
  #geom_text(vjust=1.5,alpha=0.5,size=3)+
  geom_smooth(method=lm, se=F)+
  stat_cor()+
  xlab("Module membership") + ylab("Gene significance")+
  labs(subtitle=paste("Turquoise module"))
p2<-ggplot(MMGS_brown, aes(x=x, y=y,label=label)) + 
  geom_point(shape="circle",colour="black",size=2,alpha=0.5) +
  geom_point(data=MMGS_brown %>% filter(label %in% upgenes$genes),aes(x=x,y=y),color="blue") +
  #geom_text(vjust=1.5,alpha=0.5,size=3)+
  geom_smooth(method=lm, se=F)+
  stat_cor()+
  xlab("Module membership") + ylab("Gene significance")+
  labs(subtitle=paste("Brown module"))
p3<-ggplot(MMGS_yellow, aes(x=x, y=y,label=label)) + 
  geom_point(shape="circle",colour="black",size=2,alpha=0.5) +
  geom_point(data=MMGS_yellow %>% filter(label %in% upgenes$genes),aes(x=x,y=y),color="blue") +
  #geom_text(vjust=1.5,alpha=0.5,size=3)+
  geom_smooth(method=lm, se=F)+
  stat_cor()+
  xlab("Module membership") + ylab("Gene significance")+
  labs(subtitle=paste("Yellow module"))+

layout <- "
AABB
#CC#
"ddd
p3+p2+p1 + plot_layout(design=layout) + plot_annotation(
  title = "Module Membership vs. Gene Significance",
)

# Functional Enrichment ---------------------------------------------------
library(org.Hs.eg.db)

yellow_go <- enrichGO(gene     = yellowGenes,
               OrgDb    = org.Hs.eg.db,
               keyType = "SYMBOL",
               ont      = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05,
               readable = TRUE)
turq_go <- enrichGO(gene     = turqGenes,
                      OrgDb    = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont      = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable = TRUE,
                    )
brown_go <- enrichGO(gene     = brownGenes,
                    OrgDb    = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont      = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff=0.05,
                    readable = TRUE,
                   )
allGenes<-list(turq=turqGenes,yellow=yellowGenes,brown=brownGenes)
compr<-compareCluster(allGenes,OrgDb    = org.Hs.eg.db,fun = "enrichGO",
                      keyType = "SYMBOL",
                      ont      = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff=0.05,
                      readable = TRUE)
# par(mfrow= c(1,3))
# barplot(simplify(brown_go), title="Brown module - Gene Ontology",order=T)
# barplot(simplify(yellow_go), title="Yellow module - Gene Ontology",order=T)
# barplot(simplify(turq_go), title="Turquoise module - Gene Ontology",order=T, showCategory = 6,font.size = 12)
turq_id<-bitr(turqGenes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
yellow_id<-bitr(yellowGenes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
brown_id<-bitr(brownGenes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
whole_id<-bitr(names(datExpr),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

turq_kegg<-enrichMKEGG(gene=turq_id$ENTREZID,universe = whole_id$ENTREZID)
yellow_kegg<-enrichMKEGG(gene=yellow_id$ENTREZID,universe = whole_id$ENTREZID)
brown_kegg<-enrichKEGG(gene=brown_id$ENTREZID,universe = whole_id$ENTREZID)

##EXPORT GENE LISTS---------------------------------
cool_modules<-c("yellow","turquoise","brown")
for (modules in cool_modules){
  moduleGenes<-mergedColors==modules
  moduleGenes<-names(datExpr)[moduleGenes]
  fname<-paste(modules,"_geneList.txt")
  write.table(moduleGenes,file=fname,sep="",row.names = F,col.names=F, quote=F)
}

##EXPORT NODES AND EDGES TO CYTOSCAPE---------------
modules = c("brown","turquoise","yellow");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule]
TOM.mat<-as.matrix(TOM)
modTOM = TOM.mat[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

modDEG = tibble(label=modProbes,
                DE=modProbes %>% case_match(downgenes$genes ~ "down", upgenes$genes ~ "up"),
                color=mergedColors[inModule],
                hub=modProbes %>% case_match(AllHubs~TRUE,.default = FALSE),
                hub_MM=modProbes %>% case_match(AllHubs_MM~TRUE,.default=FALSE))
modDEG<-inner_join(modDEG,bigMMlist,by="label")
##easier to do this in excel
exportNetworkToCytoscape(modTOM, edgeFile="edge_all2.txt", nodeFile = "node_all2.txt", threshold=0,nodeNames=modProbes, nodeAttr = modDEG)








#Exploratory DESeq RNA-Seq data analysis

#Load module
library(DESeq2) 
library(ggplot2)
library(pheatmap)

# Create counts_ordered from FeatureCounts
counts_ordered <- FeatureCounts[, rownames(SampleTable)]

#Expolratory Analysis accounting first for the condition parameter
dds_condition <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_ordered),
  colData = SampleTable,
  design = ~ Condition
)

dds_condition <- DESeq(dds_condition)


#Quality check by removing the dependence of the variance
vsd_condition <- vst(dds_condition, blind=TRUE) #vsd() is a variance stabilizing transformation, turns raw counts into something more like log‑expression values

#PCA plot for condtions
plotPCA(vsd, intgroup="Condition") +
  ggtitle("PCA Plot - Design: ~ Condition") +
  theme_minimal()

#Analysis accounting for mouse.genotype parameter
dds_genotype <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_ordered),
  colData = SampleTable,
  design = ~ Mouse.genotype
)

dds_genotype <- DESeq(dds_genotype)
vsd_genotype <- vst(dds_genotype, blind=TRUE)

plotPCA(vsd, intgroup="Mouse.genotype") +
  ggtitle("PCA Plot - Design: ~ Mouse.genotype") +
  theme_minimal()

#Analysis of the interaction between condition and mouse.genotype
SampleTable$group <- interaction(SampleTable$Condition,
                                 SampleTable$Mouse.genotype)#combining condition and mouse.genotype in factor "group"

dds_interaction <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_ordered),
  colData = SampleTable,
  design = ~ group
)
dds_interaction <- DESeq(dds_interaction)
vsd_interaction <- vst(dds_interaction, blind = TRUE)

plotPCA(vsd, intgroup = "group") +
  ggtitle("PCA: Condition + Genotype") +
  theme_minimal()

#Additional exploratory analysis to look into most variable genes
# select top 500 most variable genes
varGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)

mat  <- assay(vsd_interaction)[varGenes, ]
mat  <- mat - rowMeans(mat)#center genes

anno <- as.data.frame(colData(vsd_interaction)[, c("Condition", "Mouse.genotype")])
pheatmap(mat,
         annotation_col = anno,
         show_rownames = FALSE,
         main = "Top 500 most variable genes")#

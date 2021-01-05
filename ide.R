
# IDE of RNAseq data

# Inputs:
#  1. raw counts table (with correct sample names)
#  2. metadata excel file as coldata

# RNAseq_v2020


#############
# Libraries #
#############

library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(stringr)
library(readxl)


###################
# Input Variables #
###################

# path to raw counts table
raw_counts_tbl <- "/home/groups/hoolock2/u0/bd/Projects/ECP56/data/ide/raw_counts.txt"
# path to metadata excel file
# currently hardcoded to column "Group" specifying group membership for comparisons
metadata <- "/home/groups/hoolock2/u0/bd/Projects/ECP56/metadata/ECP56.notes.xlsx"
# plot file type (usually png or pdf)
plot_type <- "png"
# path to gene info file (often biomart info)
#  max counts section depends on this
#  currently hardcoded to expect 6 columns with geneID, gene.name as 1,2, and gene.description as 5
gene_info <- "/home/groups/hoolock2/u0/genomes/ensembl/mus_musculus/annotation/Mus_musculus.GRCm38.biomart_info.txt"


####################
# READ IN RAW DATA #
####################

# read in raw counts table
raw_counts <- read.table(raw_counts_tbl, header=T)

# Set the rownames of the raw_counts table to be the Gene_ID column and remove the Gene_ID column
rownames(raw_counts) <- raw_counts[,1]
raw_counts[,1] <- NULL


###########
# COLDATA #
###########

# read in coldata excel file
coldata <- as.data.frame(read_excel("/home/groups/hoolock2/u0/bd/Projects/ECP56/metadata/ECP56.notes.xlsx"))
colnames(coldata)[1] <- "sample"


####################################
# Get Total Counts for each sample #
####################################

tcounts <- colSums(raw_counts)

# make a df of total counts per sample and then plot (refer to coldata for variables)
tcounts.df <- data.frame(sample=colnames(raw_counts))
tcounts.df$total_counts <- tcounts

# export total counts df
write.table(tcounts.df, "total_counts_table.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# create barplot of total counts
myplot <- ggplot(tcounts.df, aes(x=sample,y=total_counts, fill=sample)) +
    labs(x="Sample", y="Total Gene Counts", title="Total Gene Counts per Sample after Alignment") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_bar(stat="identity") +
    theme_minimal()+
    #scale_fill_brewer(palette="Spectral") +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename=paste("total_counts_barplot", plot_type, sep="."), plot=myplot)


################################
# GET TOP 20 GENES (AVG COUNT) #
################################

# Add an average column to get the average counts for each gene
raw_counts$avg <- rowMeans(raw_counts)
# Re-add GeneID column
raw_counts$GeneID <- rownames(raw_counts)
# re-order to place GeneID first
raw_counts <- raw_counts[ ,c(ncol(raw_counts),1:ncol(raw_counts)-1)]

# get the 20 most abundant genes
max_counts <- raw_counts[order(raw_counts$avg, decreasing=T)[1:20], ]

# read in reference table that matches geneID to gene name
ref <- read.delim(gene_info, header=TRUE)
# change colnames
colnames(ref) <- c("GeneID","gene.name","source.of.gene.name","gene.type","gene.description","ensembl.family.description")
# remove any duplicate GeneID rows
ref <- ref[!duplicated(ref$GeneID), ]
# keep only GeneID, gene.name, gene.description
ref <- ref[ ,c(1,2,5)]

# merge max_counts table with subsetted ref table on "GeneID"
max_counts_final <- merge(max_counts, ref, by="GeneID")
# re-order columns
max_counts_final <- max_counts_final[ ,c(1,(ncol(max_counts_final) - 1), ncol(max_counts_final), (ncol(max_counts_final)-2), 2:(ncol(max_counts_final) -3))]
# reorder rows by avg count
max_counts_final <- max_counts_final[order(max_counts_final$avg, decreasing=T), ]

# export max counts table
write.table(max_counts_final, "max_counts.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# specify the order of the gene.name factor levels for correct x-axis order on plot
max_counts_final$gene.name <- factor(max_counts_final$gene.name, levels=unique(max_counts_final$gene.name))

# create barplot of 20 Genes with Max Counts on Avg
myplot <- ggplot(max_counts_final, aes(x=gene.name,y=avg)) +
    labs(x="Gene Name", y="Avg Gene Count", title="Genes with Highest Average Raw Counts Across All Samples") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_bar(stat="identity") +
    theme_minimal()+
    #scale_fill_brewer(palette="Spectral") +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename=paste("max_counts_barplot", plot_type, sep="."), plot=myplot)

# fix raw_counts to include GeneID as rowname and only samples as columns
raw_counts$GeneID <- NULL
raw_counts$avg <- NULL


#######################
# LOW COUNT FILTERING #
#######################

# Strategy: utilize edgeR's function "filterByExpr"

# Create the DGEList object from raw counts and group specification
y <- DGEList(counts=raw_counts, group=coldata$Group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
counts.keep <- as.data.frame(y$counts)
print(paste0("Genes remaining after low-count filtering: ", nrow(counts.keep)))


##################################
# Variance Stabilizing Transform #
##################################

# Create a dds object via DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts.keep, colData=coldata, design= ~ 0 + Group)
dds <- estimateSizeFactors(dds)

# Transform to stabilize variance
if (ncol(counts.keep) >= 30) {
  print("Using VST Transform due to sample size >= 30")
  vsd <- vst(dds, blind=FALSE)
} else {
  print("Using rlog transform due to sample size < 30")
  vsd <- rlog(dds, blind=FALSE)
}


#######################
# Sample Dist Heatmap #
#######################

# use "dist" function to calculate Euclidian distance between samples
sampleDists <- dist(t(assay(vsd)))
# visualize the distances in a heatmap
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

if (plot_type=="png") {
  png("eucl_dist_heatmap.png")
  pheatmap(sampleDistMatrix,
       clustering_distance_rows = sampleDists,
       clustering_distance_cols = sampleDists,
       col = colors)
  dev.off()
} else if (plot_type=="pdf") {
  pdf("eucl_dist_heatmap.pdf")
  pheatmap(sampleDistMatrix,
       clustering_distance_rows = sampleDists,
       clustering_distance_cols = sampleDists,
       col = colors)
  dev.off()
}


#######
# PCA #
#######

# we first select only the 500 genes showing the highest variance
#ntop = 500
#Pvars <- rowVars(assay(vsd))
#select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]

# use all genes passing low count filter for PCA
dataVST <- plotPCA(vsd, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(dataVST, "percentVar"))

# Make vst PCA Plot
myplot <- ggplot(dataVST, aes(PC1, PC2, color=Group)) +
       geom_point(size=3) +
       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
       ylab(paste0("PC2: ",percentVar[2],"% variance")) +
       ggtitle("Principal Component Analysis") +
       geom_text(aes(label=paste(name, sep="")),hjust=1, vjust=-.5) +
       theme_minimal()
ggsave(filename=paste("pca_plot.PC1PC2", plot_type, sep="."), plot=myplot)

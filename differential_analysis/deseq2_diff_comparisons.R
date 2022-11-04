#!/usr/bin/Rscript

# Use DESeq2 to perform differential expression analysis on RNAseq samples

library(edgeR)
library(limma)
library(DESeq2)

####################
# READ IN RAW DATA #
####################

# read in raw counts table
raw_counts <- read.table("/u1/bd/Projects/ECP9/ECP9.4_Sept2018_RNAseq/counts/final_counts_table.txt", sep="\t", header=TRUE)

# Set the rownames of the raw_counts table to be the Gene_ID column and remove the Gene_ID column
rownames(raw_counts) <- raw_counts[,1]
raw_counts[,1] <- NULL


##################
# Create COLDATA #
##################

# make coldata data frame to hold sample names as rows and the variables as columns (for use with specifying design for dds object)
coldata <- data.frame(name=c(colnames(raw_counts)), gr=c(rep("4", 3), rep("5", 3), rep("3", 3), rep("6", 3), rep("1", 3), rep("2", 3)),
	   					    batch=c(rep("1", 13), "2", "1", "1", "2", "1")
)

rownames(coldata) <- coldata[ ,1]


#######################
# LOW COUNT FILTERING #
#######################

countdata <- raw_counts

# Obtain CPMs
myCPM <- cpm(countdata)

# choose to keep genes with cpm > 1
thresh <- myCPM > 1

# we would like to keep genes that have at least 3 TRUES in each row of thresh
# KEEP: CPM > 1 in at least 3 samples: 14,922 genes
keep <- rowSums(thresh) >= 3
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]


################################################################################
# DIFFERENTIAL ANALYSIS #
#########################

# DESeqDataSet object #

# reference the filtered counts.keep table and coldata to make a DESeq Data Set object, define design
dds <- DESeqDataSetFromMatrix(countData=counts.keep, colData=coldata, design= ~ 0 + gr + batch)

# Estimate size factors (normalizes for library size)
dds <- estimateSizeFactors(dds)

# Run DESeq on dds object, containing filtered counts and specified design
# Design: ~ condition (no covariates or interactions)
dds <- DESeq(dds)

# to see which effects were fit with our model
resultsNames(dds)

# dispersion plot
png("dispersion_plot.png")
plotDispEsts(dds)
dev.off()

# Get results (default is last variable in design formula)
# A Wald Test (default) should be used when comparing two levels in a factor
# use the LRT test when the factor has 3+ Levels and you want to test all at once, like an ANOVA
res.2v1 <- results(dds, contrast=c("gr", "2", "1"))
res.3v1 <- results(dds, contrast=c("gr", "3", "1"))
res.3v2 <- results(dds, contrast=c("gr", "3", "2"))
res.5v4 <- results(dds, contrast=c("gr", "5", "4"))
res.6v4 <- results(dds, contrast=c("gr", "6", "4"))
res.6v5 <- results(dds, contrast=c("gr", "6", "5"))

res <- res.2v1
res <- res.3v1
res <- res.3v2
res <- res.5v4
res <- res.6v4
res <- res.6v5

mcols(res, use.names = TRUE)
summary(res)


#########
# PLOTS #
#########

# Independent Filtering
# We can observe how the number of rejections changes for various cutoffs based on mean normalized count
png("6v5.ind_filtering_plot.png")
plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of baseMean", ylab="number of rejections")
dev.off()

# p-value distribution histograms
png("6v5.pval_hist.png")
hist(res$pvalue, breaks=20, col="grey")
dev.off()

# DESeq2 MA Plot
png("6v5.MA_plot.png")
plotMA(res, main="DESeq2 MA Plot: NS_MP v FcE.Pos")
dev.off()

#############################################################################################################################

############
# ANNOTATE #
############

# get normalized counts
normcounts <- as.data.frame(counts(dds, normalized=TRUE))
normcounts$GeneID <- rownames(normcounts)

# read in Ensembl bioMart table for gene info
ref <- read.delim("/u1/genomes/Epigenetics_Core/mouse/m38/89/m38.89.biomart.info.txt", header=TRUE)
# remove duplicate rows based on GeneID
ref <- ref[!duplicated(ref[,1]), ]
colnames(ref) <- c("GeneID", "gene.name", "source.of.gene.name","gene.type","gene.description","ensembl.family.description")

# convert deseq2 results to df, set GeneID column
res <- as.data.frame(res)
res$GeneID <- rownames(res)

# merge res and normcounts to add normcounts to results
foo <- merge(res, normcounts, by="GeneID")

# merge this new df (foo) with the biomart df to add gene info
foo2 <- merge(foo, ref, by="GeneID", all.x=TRUE)

# remove unnecessary sample normcounts columns
foo2 <- foo2[ ,c(1,2,3,4,5,6,7,11,12,13,17,18,19,26,27,28,29,30)]

# export ALL table
write.table(foo2, "6v5.results.all.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# subset to sig only (padj < 0.05 and FC > 2)
foo3 <- foo2[!is.na(foo2$padj), ]
foo4 <- foo3[foo3$padj < 0.05, ]
foo5 <- foo4[abs(foo4$log2FoldChange) > 1, ]

# export SIG table
write.table(foo5, "6v5.results.sig.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)




# conda activate RNAseq_v2020

#------------------------------------------------------#

# LOAD LIBRARIES

library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(stringr)
library(readxl)

#------------------------------------------------------#

# LOAD CUSTOM FUNCTIONS

source("/home/groups/hoolock2/u0/bd/scripts_2020/RNAseq/initial_data_exploration/bulkRNA_ide_source_functions.R")

#------------------------------------------------------#

# USER VARIABLES

raw_counts_tbl <- "/home/groups/hoolock2/u0/bd/Projects/mariam_r21_mouse_ad_2022/rnaseq/data/counts/raw_counts.txt"
my_metadata <- "/home/groups/hoolock2/u0/bd/Projects/mariam_r21_mouse_ad_2022/rnaseq/metadata/r21_metadata.xlsx"
gene_info <- ""
plot_type <- "pdf"
outfile_prefix <- "/home/groups/hoolock2/u0/bd/Projects/mariam_r21_mouse_ad_2022/rnaseq/data/ide/" # run from scripts dir, set relative path for outfiles

#------------------------------------------------------#

# READ IN VARIABLE OBJECTS

# read in raw counts table
raw_counts <- read.table(raw_counts_tbl, header=T)
# Set the rownames of the raw_counts table to be the Gene_ID column and remove the Gene_ID column
rownames(raw_counts) <- raw_counts[,1]
raw_counts[,1] <- NULL

# read in metadata excel file
coldata <- as.data.frame(read_xlsx(my_metadata))

# check if colnames of raw counts table matches coldata$Sample_Name
identical(colnames(raw_counts), coldata$Sample_Name)

# read in reference table that matches geneID to gene name
ref <- read.delim(gene_info, header=TRUE)
# change colnames
colnames(ref) <- c("GeneID","gene.name","gene.type","gene.description")
# remove any duplicate GeneID rows
ref <- ref[!duplicated(ref$GeneID), ]
# keep only GeneID, gene.name, gene.description
#ref <- ref[ ,c(1,2,4)]

# create outdir if doesn't exist
dir.create(outfile_prefix)

#------------------------------------------------------#

# Plot Total Raw Counts Barplot
myplot <- plot_total_counts(raw_counts)
ggsave(filename=paste0(outfile_prefix, "total_counts_barplot.", plot_type), plot=myplot)


# Plot top 20 most abundant genes based on raw counts
myplot <- plot_top_count_genes(raw_counts, n=20, gene_info=ref)
ggsave(filename=paste0(outfile_prefix, "max_counts_barplot.", plot_type), plot=myplot)

#------------------------------------------------------#

# LOW-COUNT FILTERING

# Strategy: utilize edgeR's function "filterByExpr"

counts.keep <- lcf_edger(raw_counts, group=coldata$Group)
# export filtered counts
counts.keep.export <- counts.keep
counts.keep.export$GeneID <- rownames(counts.keep.export)
counts.keep.export <- counts.keep.export[ ,c(ncol(counts.keep.export),1:ncol(counts.keep.export)-1)]
write.table(counts.keep.export, "filtered_counts.txt", sep="\t", col.names=T, row.names=F, quote=F)

#------------------------------------------------------#

# PERFORM VARIANCE STABILIZING TRANSFORMATION

# Options of vst, rlog, or logcpm
# for now, do vst
vsd <- var_stable_transform(counts = counts.keep,
  coldata = coldata,
  method = "vst",
  formula_vars = c(0, "Group"),
  blind = FALSE)

#------------------------------------------------------#

# SAMPLE CORRELATION HEATMAP

sample_dist_heatmap(norm_counts = t(assay(vsd)),
  outfile = "eucl_dist_heatmap",
  plot_type = plot_type,
  plot_prefix = outfile_prefix)

#------------------------------------------------------#

# PCA

# use DESeq2 PCA function, plot myself with ggplot
myplot <- deseq2_pca(norm_counts = vsd, intgroup = "Group")
ggsave(filename=paste0(outfile_prefix, "pca_pc1_pc2.", plot_type), plot=myplot)

# Alternatively, do custom PCA
# use all features/genes
mypca <- perform_manual_pca(vsd = vsd, n=nrow(counts.keep))

# produce PCA biplots
myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC1",
  pc_b = "PC2",
  color_var = "Group",
  label_var = "Sample_Name")
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc1_pc2.", plot_type), plot=myplot)

myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC1",
  pc_b = "PC3",
  color_var = "Group",
  label_var = "Sample_Name")
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc1_pc3.", plot_type), plot=myplot)

myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC2",
  pc_b = "PC3",
  color_var = "Group",
  label_var = "Sample_Name")
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc2_pc3.", plot_type), plot=myplot)

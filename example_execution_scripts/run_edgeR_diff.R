
0;95;0c
# conda activate RNAseq_v2020

#------------------------------------------------------#

# LOAD LIBRARIES
library(writexl)

#------------------------------------------------------#

# LOAD CUSTOM FUNCTIONS

source("/home/groups/hoolock2/u0/bd/scripts_2020/RNAseq/differential_analysis/diff_analysis_edgeR_fxns.R")

#------------------------------------------------------#

# USER VARIABLES

# define paths to filtered counts file and metadata excel file
filtered_counts <- "../data/counts_tables/filtered_counts.txt"
my_metadata <- "../metadata/ECP67.metadata.xlsx"
#gene_info <- "/home/groups/hoolock2/u0/genomes/ensembl/papio_anubis/annotation/Panu_3.0.104.gene_info.txt"
gene_info <- "/home/groups/hoolock2/u0/bd/Projects/ECP67/data/annotation/gene_info.txt"
outfile_prefix <- "../data/diff/"

# read in data
counts.keep <- read.delim(filtered_counts, header=T)
rownames(counts.keep) <- counts.keep[ ,1]
counts.keep[ ,1] <- NULL

coldata <- as.data.frame(read_xlsx(my_metadata))

ref <- read.delim(gene_info, header=TRUE)
# change colnames
colnames(ref) <- c("GeneID", "gene.name")
#colnames(ref) <- c("GeneID","gene.name","gene.type","gene.description")
# remove any duplicate GeneID rows
ref <- ref[!duplicated(ref$GeneID), ]

my_design <- model.matrix(~ 0 + Group, data=coldata)

# create outdir if doesn't exist
dir.create(outfile_prefix)

#------------------------------------------------------#

# create DGEList object
y <- DGEList(counts=counts.keep, group=coldata$Group)
# Normalize (TMM normalization for RNA composition)
y <- calcNormFactors(y, method = "TMM")
# Get logCPM counts
logCPM <- as.data.frame(cpm(y, log=TRUE))
logCPM$GeneID <- rownames(logCPM)
logCPM2 <- logCPM[ ,c(ncol(logCPM),1:ncol(logCPM)-1)]
write.table(logCPM2, "filtered_logCPM_counts.txt", sep="\t", col.names=T, row.names=F, quote=F)

# Estimate Dispersion
y <- estimateDisp(y, design=my_design)
# fit the model using QL F-test and model design
fit <- glmQLFit(y, my_design)

#------------------------------------------------------#

# check the names of the fitted parameters
head(fit$coefficients)

# Use the fitted parameters to define the contrast matrix
my.contrasts <- makeContrasts(ko_v_ctrl = GroupKO - GroupCtrl,
  ctrlGW_v_ctrl = GroupCtrl_GW - GroupCtrl,
  koGW_v_ko = GroupKO_GW - GroupKO,
  levels = my_design)

#------------------------------------------------------#

# Perform all contrasts.   
# Create p-value distrbution histogram.   
# Create MA-plot.   
# Export results txt file of all genes used in differential analysis.

# run all contrasts and export plots and results
run_contrasts(my.contrasts=my.contrasts,
              export_path=outfile_prefix,
              gene_info=ref,
              fc=0,
              pc=0.05,
              pa=TRUE
)

# Make volcano plots
plot_volcano(infile = "../data/diff/mp_v_p.results.txt",
  gene_name_col = "gene.name",
  logFC_col = "logFC",
  padj_col = "padj",
  padj_cutoff = 0.05,
  logFC_cutoff = 0,
  int_genes = NULL,
  comparison_name = "micropatterned vs planar",
  outfile = "../data/diff/mp_vs_p.volcano_plot.pdf"
)

# export results to excel file
df <- read.delim("../data/diff/mp_v_p.results.txt")
write_xlsx(df, "../data/diff/ECP66.diff.results.xlsx")

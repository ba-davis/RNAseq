# Code related to RNAseq analysis

Typically, after FastQC, Trimming, and STAR alignment have been performed.

Need to make a raw counts table from the STAR output, explore data, differential analysis.

## Strategy of Bulk RNAseq Analysis
1. Utilize **FQC_Trim_align** scripts to perform via SLURM:
     - FastQC on raw files
     - Trimmomatic
     - FastQC on trimmed files
     - alignment and gene counts via STAR and ensembl

2. Utilize functions in **bulkRNA_ide_source_functions.R** to:
     - create raw counts table from STAR output
     - perform initial data exploration on the count data

3. Utilize functions to perform differential expression analysis
     - with DESeq2 or edgeR


## Input Requirements
**1. metadata excel file**
  - first column **must** be desired sample name (column **Sample_Name**)
  - **must** contain column **fastq** which is the fastq file prefix
    - this is the first part of the fastq file name
    - (the part preceding the general suffix ".star.ReadsPerGene.out.tab" of the STAR output
  - currently, uses column **Group** to define sample group membership (colors of PCA, groups to compare, etc)
  - rows (samples) are in a desired order that makes sense (with replicates together)
  - the sample columns of raw counts table will be in this same order

**2. raw counts table**
     - made with **compile_readcounts** function and passing metadat excel file to get sample names

**3. optional gene_info file**
     - provides matching gene names and gene descriptions for GeneIDs



Projects for reference:
  - ECP55 (non-ensembl genome and annotation, make gene info file optional?)
  - agarwal bulk RNAseq (clean up functions to be run within R, not command line)
  - new collection method of FastQC, Trimming, STAR stats?

     
TODO:
in max counts barplot function, how to correctly output txt file with different gene_info files - use descriptions?
in ide.R script, make sure export of filtered counts is in good format (include GeneID as first column?)
export cpm and logcpm counts in ide.R script? Currently exported in differential analysis

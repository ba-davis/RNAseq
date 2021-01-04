Code related to RNAseq analysis

Strategy:
Utilize FQC_Trim_align scripts to perform via SLURM:
  - FastQC on raw files
  - Trimmomatic
  - FastQC on trimmed files
  - alignment and gene counts via STAR and ensembl

Create a metadata excel file
  - first column must be desired sample name
  - must contain column "fastq" which is the fastq file prefix
  
Produce a raw counts table
  - (-d) given path to folder of ReadsPerGene.out.tab from STAR
  - (-o) desired output file name
  - (-s) number of column to use for counts (often 4 or 2)
  - (-m) path to excel file containing 1 sheet of metadata. First column must be desired sample name, also must contain a column "fastq" containing the prefix of the fastq file


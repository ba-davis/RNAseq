Code related to RNAseq analysis

Strategy:
Utilize FQC_Trim_align scripts to perform via SLURM:
  - FastQC on raw files
  - Trimmomatic
  - FastQC on trimmed files
  - alignment and gene counts via STAR and ensembl

Produce a raw counts table
  - given path to folder of ReadsPerGene.out.tab from STAR
  
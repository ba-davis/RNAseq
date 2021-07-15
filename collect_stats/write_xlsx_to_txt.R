#!/home/groups/hoolock2/u0/bd/miniconda3/envs/RNAseq_v2020/bin/Rscript

library(readxl)
library(optparse)

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="path to metadata excel file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (metadata excel file).n", call.=FALSE)
}

# store metadata file name
newfile <- gsub("xlsx", "txt", opt$metadata)

# read in the xlsx file
df <- as.data.frame(read_excel(opt$metadata))
colnames(df)[1] <- "Sample_Name"
# keep "Sample_Name" and "fastq" columns
df2 <- df[ ,c(colnames(df)[colnames(df) %in% c("Sample_Name", "fastq")])]
colnames(df2)[2] <- "Fastq_File"
df3 <- df2[ ,c(2,1)]
write.table(df3, newfile, sep="\t", col.names=T, row.names=F, quote=F)


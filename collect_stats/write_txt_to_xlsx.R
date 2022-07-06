#!/home/groups/hoolock2/u0/bd/miniconda3/envs/RNAseq_v2020/bin/Rscript

library(writexl)
library(optparse)

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="txt file to convert to xlsx", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input txt file).n", call.=FALSE)
}

# store outfile name
newfile <- gsub("txt", "xlsx", opt$infile)

# read in the txt file
df <- read.delim(opt$infile, header=T)

# output as xlsx file
write_xlsx(df, newfile)

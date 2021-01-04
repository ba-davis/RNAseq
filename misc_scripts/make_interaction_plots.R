
library("ggplot2")
library(edgeR)
library(limma)
library(DESeq2)

# read in the table of sig interaction genes
int.up <- read.delim("../synergy.table.UP.padj_0.1.txt", header=T)
int.down <- read.delim("../synergy.table.DOWN.padj_0.1.txt", header=T)

# obtain the dds object with DESeq2 by using the old diff scripts
# specifically, use the new coldata and design from the new diff analysis (~ diet + treat + diet:treat), 19,063 genes
resultsNames(dds)
# [1] "Intercept"      "diet_WSD_vs_C"  "treat_T_vs_noT" "dietWSD.treatT"

# Want to plot all 11 genes which surpassed padj < 0.1 threshhold
### choose FTL (ENSMMUG00000003909) and PHB2 (ENSMMUG00000010205) to plot

# idk, ANKYRIN REPEAT DOMAIN CONTAINING: ENSMMUG00000000494
# RYR2: ENSMMUG00000001060
# VPS13C: ENSMMUG00000001362
# FTL: ENSMMUG00000003909
# CDC42BPA: ENSMMUG00000008638
# ERC1: ENSMMUG00000010933
# NIN: ENSMMUG00000014658
# KLF8: ENSMMUG00000014678
# ZNF589: ENSMMUG00000021607

# CCT8: ENSMMUG00000003023
# DDT: ENSMMUG00000004552
# PHB2: ENSMMUG00000010205

# use DESeq2 plotCounts function to obtain a normalized count value
# normalizes counts by sequencing depth and adds a pseudocount of 1/2 to allow for log scale plotting
gene="ENSMMUG00000000494"
d <- plotCounts(dds, gene=gene, intgroup=c("diet", "treat"), returnData=TRUE)

# add a new column
d$group <- paste(d$diet, d$treat, sep=".")
d$connect <- c(1:24)

# obtain the centroids of the 4 data classes (mean)
centroids <- aggregate(count~group, d, mean)

# add a diet and treat column to centroids df for plotting colors
centroids$diet <- c("C", "C", "WSD", "WSD")
centroids$treat <- c("noT", "T", "noT", "T")
centroids$connect <- c("red", "red", "blue", "blue")

my.title <- "Homo ANKRD Expression"
my.file <- "Homo_ANKRD.intPlot.png"

myplot <- ggplot(d, aes(x=treat, y=count, color=diet)) +
       geom_point(position=position_jitter(w=0.1,h=0)) +
       geom_point(data=centroids, size=5) +
       geom_line(data=centroids, aes(group=connect)) +
       xlab("Testosterone Treatment") +
       ggtitle(my.title) +
       #scale_y_log10(breaks=c(5000,10000,15000)) +
       scale_color_manual(values=c("C"="darkviolet", "WSD"="green2")) +
       theme_light()
ggsave(filename=my.file, plot=myplot)

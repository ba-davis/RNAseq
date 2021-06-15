
# Make a volcano plot from an RNAseq differential analysis

library(ggplot2)
library(ggrepel)

myfile <- "/home/groups/hoolock2/u0/bd/Projects/chavez/rhesus_fragment_solo/data/diff/try1/Frag4C_vs_Unfrag4C/Frag_4C_vs_Unfrag_4C.results.txt"

# read in the full diff analysis output
mydat <- read.delim(myfile, header=T)

#---------------------------------------------------------------------------------------------------
# subset to -log10(adjusted pvalue) and Log2FC columns and geneID and gene name columns
mydat2 <- mydat[ ,c(2,6,1,39)]
mydat2$padj <- -log10(mydat2$padj)

# add a colors column to give colors to points based on -log10pval and logFC direction
# variables to allow specific thresholds:
#   padj: desired padj sig cutoff (default 0.05)
#   logFC: desired log2FC value for sig cutoff (default -1 and 1)
# if -log10(pval) < -log10(padj), color the points gray
# if -1og10(pval) > -log10(padj) AND logFC < -(logFC), color the points blue
# if -1og10(pval) > -log10(padj) AND logFC > logFC, color the points red

my_padj <- 0.05
my_logFC <- 1
my_genes <- c("NSL1", "BRK1","ESCO1", "OFD1", "WEE2")

mydat2$key <- "not_sig"
mydat2[which(mydat2$padj > -log10(my_padj) & mydat2$logFC < -(my_logFC)),"key"] <- "sig_down"
mydat2[which(mydat2$padj > -log10(my_padj) & mydat2$logFC > my_logFC),"key"] <- "sig_up"


#-----------------------------------------------------------------------------------------------------
# Volcano Plot
myplot <- ggplot(mydat2, aes(logFC, padj, group=key)) +
  geom_point(aes(color=key), size=1.5) +
  xlab("logFC") +
  ylab("-log10(P-val)") +
  ggtitle("Frag4C vs Unfrag4C Volcano Plot") +
  #geom_text(aes(label=ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),'')),hjust=0,vjust=0) +
  #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  geom_label_repel(aes(label = ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),''),
                  fill = key),
		  color = "white",
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50',
		  min.segment.length = 0) +
  scale_color_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
  scale_fill_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
  guides(fill=FALSE) +
  theme_classic()
ggsave(filename="volcano_plot.pdf")


# no box around labels plot
myplot <- ggplot(mydat2, aes(logFC, padj, group=key)) +
  geom_point(aes(color=key), size=1.5) +
  xlab("logFC") +
  ylab("-log10(P-val)") +
  ggtitle("Frag4C vs Unfrag4C Volcano Plot") +
  #geom_text(aes(label=ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),'')),hjust=0,vjust=0) +
  #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  #geom_label_repel(aes(label = ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),''),
  #                fill = key),
  #                color = "white",
  #                box.padding   = 0.35,
  #                point.padding = 0.5,
  #                segment.color = 'grey50',
  #                min.segment.length = 0) +
  geom_text_repel(aes(label = ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),'')),
                  min.segment.length = 0) +
  scale_color_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
  #scale_fill_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
  #guides(fill=FALSE) +
  theme_classic()
ggsave(filename="volcano_plot.noBox..pdf")
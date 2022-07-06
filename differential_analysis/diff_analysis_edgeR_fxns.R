

# edgeR for differential analysis functions

library(edgeR)
library(readxl)
library(ggplot2)
library(ggrepel)

# set max overlaps for ggrepel to infinite for volcano plot gene labeling
options(ggrepel.max.overlaps = Inf)

#-----------------------------------------------------------------------------------------#
# FUNCTIONS
# 1. edgeStats()
# 2. edgeMAplot()
# 3. three_lines()
# 4. run_contrasts()
# 5. plot_volcano()

#-----------------------------------------------------------------------------------------#


# Inputs:
# 1. Filtered counts (NOT variance stabilized, raw counts just low count filtered)
# 2. coldata data.frame (made from metadata excel file)
# 3.


#-----------------------------------------------------------------#

# edgeR sumStats function
edgeStats <- function(x,fc=1,pc=0.05,pa=FALSE){
  # Use the adjusted p-value for the final two columns and the upfc,dnfc cutoff
  up <- x$table$logFC > 0
  dn <- x$table$logFC < 0

  padj <- p.adjust(x$table$PValue, method="BH") < pc
  pval <- x$table$PValue < pc
  if(pa){
    ng.up <- x$table$logFC > fc & padj
    ng.dn <- x$table$logFC < -1*fc & padj
  }else{
    ng.up <- x$table$logFC[x$table$logFC > fc & pval]
    ng.dn <- x$table$logFC[x$table$logFC < -1*fc & pval]
  }
  rmin <- min(x$table$logFC[pval & up], na.rm=TRUE)
  rmax <- max(x$table$logFC[pval & dn], na.rm=TRUE)

  cat(c("padj", "pval","min.up","min.dn", "up","dn"),"\n")
  cat(c(sum(padj,na.rm=TRUE), sum(pval,na.rm=TRUE), sprintf("%.3f",rmin), sprintf("%.3f",rmax), sum(ng.up, na.rm=TRUE), sum(ng.dn, na.rm=TRUE)), "\n")
  return(list(up=which(pval & up), dn=which(pval & dn), upfc=which(ng.up), dnfc=which(ng.dn), padj=which(padj)))
}

# edgeR MA plot function
edgeMAplot <- function(x, goi.list, y.ax="log2FC", ptitle="MA plot", ymin=-10, ymax=10){
  plot(x$table$logCPM, x$table$logFC, ylim=c(ymin,ymax), xlab="mean of logCPM", ylab=y.ax, main=ptitle, pch=20, col="dimgrey", cex=0.4)
  abline(h=c(-1,1), col="grey", lty="dashed")
  points(x$table$logCPM[goi.list$upfc], x$table$logFC[goi.list$upfc], col="firebrick2", pch=20, cex=0.9)
  points(x$table$logCPM[goi.list$dnfc], x$table$logFC[goi.list$dnfc], col="royalblue3", pch=20, cex=0.9)
  legend("bottomright", legend=c(paste("up (",length(goi.list$upfc),")",sep=""),
                                 paste("down (",length(goi.list$dnfc),")",sep="")), col=c("firebrick2","royalblue3"), pch=19)
  three_lines(cols=c("black","black","black"))
}
# 3 lines function
three_lines <- function(hp=c(1,0,-1),style=c(2,1,2),cols=c("grey","grey","grey")){
  abline(h=hp[1], col=cols[1], lty=style[1])
  abline(h=hp[2], col=cols[2], lty=style[2])
  abline(h=hp[3], col=cols[3], lty=style[3])
}

#------------------------------------------------------------------------------------------------------------------#

# function to run DE contrast and export results
# expects a "my.contrasts" matrix of desired contrasts
# inputs:
#   my.contrasts: edgeR my.contrasts matrix
#   export_path: prefix to place exported files
#   gene_info: cleaned gene_info data.frame, first column is GeneID
#   norm_counts: include normalized CPM counts in output (default TRUE)
#   log: logCPM the normalized counts (default TRUE)
#   fc: L2FC cutoff for sig genes in MA plot (default 1 = Fold Change of 2)
#   pc: pvalue cutoff for sig genes in MA plot (default 0.05)
#   pa: whether to use adjusted pval for sig genes in MA plot (default TRUE)
run_contrasts <- function(my.contrasts, export_path="./", gene_info, norm_counts=TRUE, log=TRUE, fc=1, pc=0.05, pa=TRUE) {

  # loop through my.contrasts
  for (i in 1:length(colnames(my.contrasts))) {
    # grab name of contrast
    resname <- colnames(my.contrasts)[i]
    print(paste0("Performing contrast ", colnames(my.contrasts)[i]))
    # run the QLF test
    res <- glmQLFTest(fit, contrast=my.contrasts[,colnames(my.contrasts)[i]])

    # make pval distribution plot
    png(paste0(export_path, resname, ".pvalue.dist.png"))
    hist(res$table$PValue)
    dev.off()
    hist(res$table$PValue)

    # MA Plot
    my.stats <- edgeStats(res, fc=fc, pc=pc, pa=pa)
    png(paste0(export_path, resname, ".my_ma_plot.png"))
    edgeMAplot(res, my.stats, y.ax="log2FC", ptitle=paste0("MA Plot:\n", resname, "\n", "Sig Genes = L2FC > ", fc, " and padj < ", pc), ymin=-10, ymax=10)
    dev.off()
    edgeMAplot(res, my.stats, y.ax="log2FC", ptitle=paste0("MA Plot:\n", resname, "\n", "Sig Genes = L2FC > ", fc, " and padj < ", pc), ymin=-10, ymax=10)

    # Make export table
    mytab <- res$table
    mytab$padj <- p.adjust(mytab$PValue, method="BH")
    mytab$GeneID <- rownames(mytab)
    # merge DE results with normalized counts
    #if (log==TRUE) {
    #  mydat <- merge(mytab, tmm.log, by="GeneID")
    #}
    #if (log==FALSE) {
    #  mydat <- merge(mytab, tmm, by="GeneID")
    #}

    # merge again with annotation info
    myres <- merge(mytab, gene_info, by="GeneID")

    # export diff results
    write.table(myres, paste0(export_path, resname, ".results.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
}

#-----------------------------------------------------------------------------------------------------------------------------------#

# INPUTS:        infile: text file with columns for: gene.name, logFC, padj
#         gene_name_col: name of column in infile holding gene names
#             logFC_col: name of column in infile holding log fold changes
#              padj_col: name of column in infile holding adjusted p-values
#           padj_cutoff: value to use for padj cutoff to define sig DE genes
#          logFC_cutoff: value to use for logFC cutoff to define sig DE genes
#             int_genes: genes of interest to label
#       comparison_name: name of comparison for title of plot
#               outfile: desired name of plot file
plot_volcano <- function(infile, gene_name_col="gene.name", logFC_col="logFC", padj_col="padj", padj_cutoff=0.05, logFC_cutoff=1, int_genes=NULL, comparison_name, outfile="volcano_plot.pdf") {
  # read in infile
  mydat <- read.delim(infile, header=T)

  # change necessary colnames
  colnames(mydat)[colnames(mydat)==gene_name_col] <- "gene.name"
  colnames(mydat)[colnames(mydat)==logFC_col] <- "logFC"
  colnames(mydat)[colnames(mydat)==padj_col] <- "padj"

  # Take -log10 of the adjusted p-value
  mydat$padj <- -log10(mydat$padj)

  # Create new column for plot colors based on significance
  mydat$key <- "not_sig"
  mydat[which(mydat$padj > -log10(padj_cutoff) & mydat$logFC < -(logFC_cutoff)),"key"] <- "sig_down"
  mydat[which(mydat$padj > -log10(padj_cutoff) & mydat$logFC > logFC_cutoff),"key"] <- "sig_up"

  # Plot with ggplot
  p <- ggplot(mydat, aes(logFC, padj, group=key)) +
    geom_point(aes(color=key), size=1.5) +
    xlab("logFC") +
    ylab("-log10(P-val)") +
    ggtitle(paste0(comparison_name, " Volcano Plot")) +
    #geom_text(aes(label=ifelse(external_gene_name %in% my_genes, as.character(external_gene_name),'')),hjust=0,vjust=0) +
    #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
    geom_label_repel(aes(label = ifelse(gene.name %in% int_genes, as.character(gene.name),''),
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
  ggsave(filename=outfile)

  return(p)
}

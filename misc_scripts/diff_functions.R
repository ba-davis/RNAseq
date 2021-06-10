

# function to run DE contrast and export results
# expects a "my.contrasts" matrix of desired contrasts
# inputs:
#   my.contrasts: edgeR my.contrasts matrix
#   export_path: prefix to place exported files
#   gene_info: biomart info, first columns is GeneIDs (header=T)
#   norm_counts: include normalized CPM counts in output (default TRUE)
#   log: logCPM the normalized counts (default TRUE)
#   fc: L2FC cutoff for sig genes in MA plot (default 1 = Fold Change of 2)
#   pc: pvalue cutoff for sig genes in MA plot (default 0.05)
#   pa: whether to use adjusted pval for sig genes in MA plot (default TRUE)
run_contrasts <- function(my.contrasts, export_path="../data/diff/", gene_info, norm_counts=TRUE, log=TRUE, fc=1, pc=0.05, pa=TRUE) {
  # read in the gene info file
  annot <- read.delim(gene_info, header=T)
  colnames(annot)[1] <- "GeneID"
  annot <- annot[!duplicated(annot$GeneID), ]

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

    # MA Plot
    #source('/home/groups/hoolock2/u0/bd/my_scripts/RNAseq/DE_ref_funcs.R')
    my.stats <- edgeStats(res, fc=fc, pc=pc, pa=pa)
    png(paste0(export_path, resname, ".my_ma_plot.png"))
    edgeMAplot(res, my.stats, y.ax="log2FC", ptitle=paste0("MA Plot:\n", resname, "\n", "Sig Gene Cutoffs: L2FC: ", fc, " and padj < ", pc), ymin=-10, ymax=10)
    dev.off()

    # Make export table
    mytab <- res$table
    mytab$padj <- p.adjust(mytab$PValue, method="BH")
    mytab$GeneID <- rownames(mytab)

    print(paste0("There are ", nrow(mytab[mytab$padj < pc, ]), " genes with padj < ", pc))
    
    # merge DE results with normalized counts
    if (log==TRUE) {
      mydat <- merge(mytab, tmm.log, by="GeneID")
    }
    if (log==FALSE) {
      mydat <- merge(mytab, tmm, by="GeneID")
    }

    # merge again with annotation info
    myres <- merge(mydat, annot, by="GeneID")

    # export diff results
    write.table(myres, paste0(export_path, resname, ".results.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

    # get sig results
    myres.sig <- myres[myres$padj < pc, ]
    # export sig diff results
    write.table(myres.sig, paste0(export_path, resname, ".sig_results.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
}

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
# 3 lines function
three_lines <- function(hp=c(1,0,-1),style=c(2,1,2),cols=c("grey","grey","grey")){
  abline(h=hp[1], col=cols[1], lty=style[1])
  abline(h=hp[2], col=cols[2], lty=style[2])
  abline(h=hp[3], col=cols[3], lty=style[3])
}

# DESeq2 sumStats function
sumStats <- function(x,fc=1,pc=0.05,pa=FALSE){
  # Use the adjusted p-value for the final two columns and the upfc,dnfc cutoff
  up <- x$log2FoldChange > 0
  dn <- x$log2FoldChange < 0

  padj <- x$padj < pc
  pval <- x$pvalue < pc
  if(pa){
    ng.up <- x$log2FoldChange > fc & padj
    ng.dn <- x$log2FoldChange < -1*fc & padj
  }else{
    ng.up <- x$log2FoldChange[x$log2FoldChange > fc & pval]
    ng.dn <- x$log2FoldChange[x$log2FoldChange < -1*fc & pval]
  }
  rmin <- min(x$log2FoldChange[pval & up], na.rm=TRUE)
  rmax <- max(x$log2FoldChange[pval & dn], na.rm=TRUE)

  cat(c("padj", "pval","min.up","min.dn", "up","dn"),"\n")
  cat(c(sum(padj,na.rm=TRUE), sum(pval,na.rm=TRUE), sprintf("%.3f",rmin), sprintf("%.3f",rmax), sum(ng.up, na.rm=TRUE), sum(ng.dn, na.rm=TRUE)), "\n")
  return(list(up=which(pval & up), dn=which(pval & dn), upfc=which(ng.up), dnfc=which(ng.dn), padj=which(padj)))
}

# DESeq2 MA Plot Function
splot <- function(x,goi.list, y.ax="log2FC", ptitle="MA plot", ymin=-10, ymax=10){
  plot(x$baseMean, x$log2FoldChange, ylim=c(ymin,ymax), log="x", xlab="mean of normalized counts", ylab=y.ax, main=ptitle, pch=20, col="dimgrey", cex=0.4)
  abline(h=c(-1,1), col="grey", lty="dashed")
  points(x$baseMean[goi.list$upfc], x$log2FoldChange[goi.list$upfc], col="firebrick2", pch=20, cex=0.9)
  points(x$baseMean[goi.list$dnfc], x$log2FoldChange[goi.list$dnfc], col="royalblue3", pch=20, cex=0.9)
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

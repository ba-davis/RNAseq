

# conda activate RNAseq_v2020

library(ggplot2)

limma_ma_plot <- function(infile, comparison_name, uplfc_line=1, dnlfc_line=-1, ylim_min=-6, ylim_max=6, file_type="pdf") {
  df <- read.delim(infile, header=T)

  df$genes <- "not_sig"
  df$genes[df$adj.P.Val < 0.05 & df$logFC > 0] <- "sig_up"
  df$genes[df$adj.P.Val < 0.05 & df$logFC < 0] <- "sig_down"

  #df$genes <- factor(df$genes, levels=c("not_sig", "sig_up", "sig_down"))
  
  p <- ggplot(df, aes(AveExpr, logFC, group=genes)) +
    geom_point(aes(color=genes), size=1.5) +
    xlab("Mean Expression") +
    ylab("log2FC") +
    ggtitle(paste0(comparison_name, " MA Plot")) +
    scale_color_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
    scale_fill_manual(values=c(not_sig="gray64", sig_up="firebrick2", sig_down="cornflowerblue")) +
    guides(fill=FALSE) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = uplfc_line, linetype="dotted") +
    geom_hline(yintercept = dnlfc_line, linetype="dotted") +
    ylim(ylim_min, ylim_max)
  ggsave(filename=paste0(comparison_name, "_ma_plot.", file_type))
}


count_sigs <- function(infile) {
  df <- read.delim(infile, header=T)

  df.sig <- df[!(is.na(df$adj.P.Val)), ]
  df.up <- df.sig[df.sig$adj.P.Val < 0.05 & df.sig$logFC > 0, ]
  df.dn <- df.sig[df.sig$adj.P.Val < 0.05 & df.sig$logFC < 0, ]

  print(infile)
  print(paste0("sig up genes: ", nrow(df.up)))
  print(paste0("sig down genes: ", nrow(df.dn)))
}

# EXAMPLE USAGE:
# limma_ma_plot("c13.diff.results.txt",
#   comparison_name="c13",
#   ylim_min=-6,
#   ylim_max=6,
#   file_type="pdf")

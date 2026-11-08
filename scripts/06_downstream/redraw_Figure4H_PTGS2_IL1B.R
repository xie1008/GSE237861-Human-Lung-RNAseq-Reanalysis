#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
script_path <- sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grep("^--file=", commandArgs(trailingOnly=FALSE))][1])
project <- normalizePath(file.path(dirname(script_path), "../.."), mustWork=TRUE)
out <- file.path(project, "results")

# Use the already-generated ROC coordinates and fixed AUC estimates; no statistics are recalculated here.
auc_labels <- c(PTGS2="PTGS2 (AUC=0.245)", IL1B="IL1B (AUC=0.684)")
files <- setNames(file.path(out, c("ROC_coordinates_PTGS2_14samples.tsv", "ROC_coordinates_IL1B_14samples.tsv")), names(auc_labels))
curves <- do.call(rbind, lapply(names(files), function(g) {
  d <- read.delim(files[[g]], check.names=FALSE)
  d <- d[is.finite(d$threshold) & is.finite(d$sensitivity) & is.finite(d$specificity),]
  d$FPR <- 1 - d$specificity
  d$gene <- g
  d <- d[order(d$FPR, d$sensitivity),]
  # Remove threshold-order backtracking and retain a monotone ROC envelope.
  d$TPR <- cummax(d$sensitivity)
  d <- d[!duplicated(d$FPR, fromLast=TRUE), c("gene","FPR","TPR")]
  d
}))
curves$gene_label <- unname(auc_labels[curves$gene])
write.table(curves, file.path(out, "ROC_plot_data_main14_monotone.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

p <- ggplot(curves, aes(FPR, TPR, color=gene_label, group=gene_label)) +
  geom_step(direction="vh", linewidth=0.9) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey70", linewidth=0.55) +
  coord_equal(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  scale_color_manual(values=c("PTGS2 (AUC=0.245)"="#0072B2", "IL1B (AUC=0.684)"="#D55E00"), name=NULL) +
  labs(x="1 - specificity", y="Sensitivity", title="PTGS2 and IL1B ROC curves") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=12), legend.position="right", legend.text=element_text(size=10), axis.title=element_text(size=11), axis.text=element_text(size=10))
ggsave(file.path(out,"Figure_4H_PTGS2_IL1B_ROC.pdf"), p, width=5.6, height=4.6, useDingbats=FALSE)
ggsave(file.path(out,"Figure_4H_PTGS2_IL1B_ROC.tiff"), p, width=5.6, height=4.6, dpi=600, compression="lzw")
cat("Redrawn Figure 4H without ROC backtracking segments. AUC values were not recalculated.\n")

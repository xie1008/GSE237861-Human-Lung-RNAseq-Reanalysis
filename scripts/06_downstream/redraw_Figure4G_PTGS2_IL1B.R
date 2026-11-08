#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
script_path <- sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grep("^--file=", commandArgs(trailingOnly=FALSE))][1])
project <- normalizePath(file.path(dirname(script_path), "../.."), mustWork=TRUE)
out <- file.path(project, "results")

# Plot-only redraw from the saved normalized expression table; no statistics are recalculated.
d <- read.delim(file.path(out, "PTGS2_IL1B_TMM_log2CPM_long.tsv"), check.names=FALSE)
d$group <- factor(d$group, levels=c("Control","Sepsis"))
d$gene <- factor(d$gene, levels=c("PTGS2","IL1B"))

p <- ggplot(d, aes(x=group, y=log2CPM, color=group)) +
  geom_boxplot(width=0.48, linewidth=0.65, outlier.shape=NA, fill="white", alpha=1) +
  geom_jitter(width=0.075, height=0, size=1.75, alpha=0.9) +
  facet_wrap(~gene, nrow=1, scales="free_y") +
  scale_color_manual(values=c(Control="#666666", Sepsis="#E66101")) +
  labs(x=NULL, y="Expression", title="PTGS2 and IL1B expression on the verified processed scale") +
  theme_classic(base_size=10.5) +
  theme(
    plot.title=element_text(size=11, face="plain", hjust=0),
    strip.background=element_rect(fill="white", color="black", linewidth=0.6),
    strip.text=element_text(size=9.5, face="plain"),
    axis.title=element_text(size=10),
    axis.text=element_text(size=9, color="black"),
    axis.line=element_line(color="black", linewidth=0.55),
    axis.ticks=element_line(color="black", linewidth=0.45),
    legend.position="none",
    panel.spacing=unit(0.75,"lines")
  )

ggsave(file.path(out,"Figure_4G_PTGS2_IL1B.pdf"), p, width=7.2, height=3.8, useDingbats=FALSE)
ggsave(file.path(out,"Figure_4G_PTGS2_IL1B.tiff"), p, width=7.2, height=3.8, dpi=600, compression="lzw")
cat("Redrawn Figure 4G using saved TMM log2-CPM values; no statistics recalculated.\n")

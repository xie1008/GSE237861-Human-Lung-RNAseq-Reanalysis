#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(pROC)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grep("^--file=", commandArgs(trailingOnly=FALSE))][1])
project <- normalizePath(file.path(dirname(script_path), "../.."), mustWork=TRUE)
out <- file.path(project, "results")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
counts_file <- file.path(project, "data/GSE237861_lung14_gene_counts.tsv")
ann_file <- file.path(project, "data/GSE237861_Ensembl112_gene_annotation.tsv")
manifest_file <- file.path(project, "metadata/lung_14sample_manifest.tsv")

counts <- read.delim(counts_file, check.names=FALSE, stringsAsFactors=FALSE)
ann <- read.delim(ann_file, check.names=FALSE, stringsAsFactors=FALSE)
manifest <- read.delim(manifest_file, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(counts)==63140, ncol(counts)==15, nrow(manifest)==14)
rownames(counts) <- counts$gene_id
count_mat <- as.matrix(counts[, manifest$GSM, drop=FALSE])
storage.mode(count_mat) <- "integer"
rownames(manifest) <- manifest$GSM
manifest$group <- factor(manifest$group, levels=c("Control","Sepsis"))

targets <- c(PTGS2="ENSG00000073756", IL1B="ENSG00000125538")
target_ann <- ann[match(unname(targets), ann$gene_id),]
if(any(is.na(target_ann$gene_id)) || any(duplicated(target_ann$gene_id)) || !all(target_ann$gene_name==names(targets))) stop("Target mapping failed")
write.table(target_ann, file.path(out,"target_gene_annotation.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

write.table(data.frame(GSM=manifest$GSM,SRX=manifest$SRX,SRR=manifest$SRR,group=manifest$group,tissue=manifest$tissue,title=manifest$title), file.path(out,"sample_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

fit_set <- function(samples, label) {
  yfull <- DGEList(counts=count_mat[,samples,drop=FALSE], group=manifest[samples,"group"])
  yfull <- calcNormFactors(yfull, method="TMM")
  full_lcpm <- cpm(yfull, log=TRUE, prior.count=2)
  keep <- filterByExpr(yfull, design=model.matrix(~yfull$samples$group))
  y <- yfull[keep,,keep.lib.sizes=FALSE]
  design <- model.matrix(~group, data=y$samples)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  qlf <- glmQLFTest(fit, coef="groupSepsis")
  tab <- topTags(qlf, n=Inf, sort.by="none")$table
  tab$gene_id <- rownames(tab); tab$gene_name <- ann$gene_name[match(tab$gene_id, ann$gene_id)]
  tab <- tab[,c("gene_id","gene_name","logFC","logCPM","F","PValue","FDR")]
  write.table(tab, file.path(out,paste0("edgeR_QL_results_",label,".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  norm <- data.frame(GSM=colnames(yfull), group=yfull$samples$group, lib.size=yfull$samples$lib.size, norm.factor=yfull$samples$norm.factors, effective.lib.size=yfull$samples$lib.size*yfull$samples$norm.factors)
  write.table(norm, file.path(out,paste0("TMM_normalization_",label,".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(data.frame(gene_id=rownames(full_lcpm), full_lcpm, check.names=FALSE), file.path(out,paste0("log2CPM_",label,".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  list(y=y, lcpm=full_lcpm, tab=tab, keep=keep)
}

main <- fit_set(colnames(count_mat), "main14")
write.table(data.frame(gene_id=rownames(count_mat), count_mat, check.names=FALSE), file.path(out,"target_and_all_raw_counts_input.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
target_counts <- data.frame(gene=names(targets), gene_id=unname(targets), count_mat[unname(targets),,drop=FALSE], check.names=FALSE)
write.table(target_counts, file.path(out,"PTGS2_IL1B_raw_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

expression <- data.frame()
for(g in names(targets)) {
  gid <- targets[[g]]
  expression <- rbind(expression, data.frame(gene=g, gene_id=gid, GSM=colnames(main$lcpm), group=manifest[colnames(main$lcpm),"group"], log2CPM=as.numeric(main$lcpm[gid,]), raw_count=count_mat[gid,]))
}
write.table(expression, file.path(out,"PTGS2_IL1B_TMM_log2CPM_long.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

mw <- do.call(rbind, lapply(names(targets), function(g) {
  z <- expression[expression$gene==g,]
  x <- z$log2CPM[z$group=="Sepsis"]; y <- z$log2CPM[z$group=="Control"]
  data.frame(gene=g, method="Mann-Whitney U", n_sepsis=length(x), n_control=length(y), median_sepsis=median(x), median_control=median(y), p_value=wilcox.test(x,y,exact=FALSE)$p.value)
}))
edge_targets <- main$tab[match(unname(targets),main$tab$gene_id),]
diff <- data.frame(gene=names(targets), gene_id=unname(targets), log2FC=edge_targets$logFC, edgeR_PValue=edge_targets$PValue, edgeR_FDR=edge_targets$FDR)
diff <- merge(diff,mw,by=c("gene"),sort=FALSE)
write.table(diff, file.path(out,"PTGS2_IL1B_differential_expression.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

roc_one <- function(expr, samples_label, label) {
  keep <- samples_label
  response <- factor(as.character(manifest[keep,"group"]), levels=c("Control","Sepsis"))
  pred <- if(is.null(names(expr))) as.numeric(expr) else as.numeric(expr[keep])
  if(length(unique(keep)) != length(keep) || any(is.na(response))) stop(paste("ROC sample mapping failed", label, paste(keep, collapse=",")))
  r <- roc(response, pred, levels=c("Control","Sepsis"), direction="<", quiet=TRUE)
  ci <- ci.auc(r, method="bootstrap", boot.n=2000, stratified=TRUE, progress="none")
  coords <- coords(r, x="all", ret=c("threshold","sensitivity","specificity","ppv","npv"), transpose=FALSE)
  coords$gene <- label; coords$direction_fixed <- "higher expression predicts Sepsis"; coords <- coords[,c("gene","direction_fixed","threshold","sensitivity","specificity","ppv","npv")]
  write.table(coords, file.path(out,paste0("ROC_coordinates_",label,"_",length(keep),"samples.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  data.frame(gene=label, n=length(keep), AUC=as.numeric(auc(r)), CI95_low=ci[1], CI95_high=ci[3], direction="higher expression predicts Sepsis")
}

roc_main <- do.call(rbind,lapply(names(targets), function(g) roc_one(expression$log2CPM[expression$gene==g], colnames(main$lcpm), g)))
write.table(roc_main, file.path(out,"ROC_AUC_main14.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

excluded <- list(main14=character(), exclude_SRR25368433="SRR25368433", exclude_SRR25368428="SRR25368428", exclude_both=c("SRR25368433","SRR25368428"))
sens_roc <- do.call(rbind,lapply(names(excluded), function(label) {
  excluded_gsm <- manifest[manifest$SRR %in% excluded[[label]],"GSM"]
  ss <- setdiff(colnames(count_mat), excluded_gsm)
  f <- fit_set(ss,label)
  do.call(rbind,lapply(names(targets), function(g) roc_one(f$lcpm[targets[[g]],], ss, paste0(g,"_",label))))
}))
write.table(sens_roc, file.path(out,"ROC_AUC_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

sens_diff <- do.call(rbind,lapply(names(excluded), function(label) {
  excluded_gsm <- manifest[manifest$SRR %in% excluded[[label]],"GSM"]
  ss <- setdiff(colnames(count_mat), excluded_gsm)
  f <- fit_set(ss,label)
  out_rows <- lapply(names(targets), function(g) {
    gid <- targets[[g]]; z <- f$lcpm[gid,]; grp <- manifest[names(z),"group"]
    mwp <- wilcox.test(z[grp=="Sepsis"],z[grp=="Control"],exact=FALSE)$p.value
    hit <- f$tab[f$tab$gene_id==gid,,drop=FALSE]
    data.frame(sensitivity_set=label,gene=g,gene_id=gid,n=length(ss),log2FC=if(nrow(hit))hit$logFC else NA_real_,edgeR_PValue=if(nrow(hit))hit$PValue else NA_real_,edgeR_FDR=if(nrow(hit))hit$FDR else NA_real_,MannWhitney_PValue=mwp,filter_retained=gid %in% rownames(f$y))
  }); do.call(rbind,out_rows)
}))
write.table(sens_diff, file.path(out,"differential_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

plot_expr <- ggplot(expression, aes(x=group,y=log2CPM, color=group, fill=group)) + geom_boxplot(width=.55, alpha=.15, outlier.shape=NA) + geom_jitter(width=.09, size=2.2, alpha=.9) + facet_wrap(~gene, scales="free_y", nrow=1) + scale_color_manual(values=c(Control="#4C78A8",Sepsis="#D62728")) + scale_fill_manual(values=c(Control="#4C78A8",Sepsis="#D62728")) + labs(x=NULL,y="TMM-normalized log2-CPM",title="GSE237861 lung external cohort") + theme_classic(base_size=12) + theme(legend.position="top")
ggsave(file.path(out,"Figure_4G_PTGS2_IL1B.pdf"),plot_expr,width=7.2,height=3.8,useDingbats=FALSE)
ggsave(file.path(out,"Figure_4G_PTGS2_IL1B.tiff"),plot_expr,width=7.2,height=3.8,dpi=600,compression="lzw")

roc_plot_data <- do.call(rbind,lapply(names(targets), function(g) { z <- expression[expression$gene==g,]; r <- roc(z$group,z$log2CPM,levels=c("Control","Sepsis"),direction="<",quiet=TRUE); data.frame(gene=g,FPR=1-r$specificities,TPR=r$sensitivities) }))
write.table(roc_plot_data,file.path(out,"ROC_plot_data_main14.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
roc_plot <- ggplot(roc_plot_data,aes(FPR,TPR,color=gene)) + geom_step(direction="vh",linewidth=.9) + geom_abline(slope=1,intercept=0,linetype=2,color="grey50") + coord_equal() + scale_color_manual(values=c(PTGS2="#1B9E77",IL1B="#D95F02")) + labs(x="1 - specificity",y="sensitivity",color=NULL,title="Apparent discrimination in GSE237861 lung samples") + theme_classic(base_size=12) + theme(legend.position="top")
ggsave(file.path(out,"Figure_4H_PTGS2_IL1B_ROC.pdf"),roc_plot,width=5.2,height=4.6,useDingbats=FALSE)
ggsave(file.path(out,"Figure_4H_PTGS2_IL1B_ROC.tiff"),roc_plot,width=5.2,height=4.6,dpi=600,compression="lzw")

writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(capture.output({cat("edgeR",as.character(packageVersion("edgeR")),"\\n");cat("pROC",as.character(packageVersion("pROC")),"\\n");cat("ggplot2",as.character(packageVersion("ggplot2")),"\\n")}),file.path(out,"package_versions.txt"))
cat("Completed",out,"\\n")

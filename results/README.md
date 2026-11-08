# GSE237861 PTGS2/IL1B external lung validation

This directory contains the final two-gene external-cohort analysis using the existing 14-sample Ensembl Release 112 featureCounts matrix. No FASTQ, trimming, alignment, or counting was repeated.

## Gene identity

- PTGS2: ENSG00000073756
- IL1B: ENSG00000125538

Both IDs are unique in the Ensembl Release 112 annotation. PTGS2 and IL1B raw counts are present for all 14 samples. IL1B has low counts and was removed by edgeR `filterByExpr` in the main 14-sample differential model; its normalized expression and ROC analysis were retained without imputation. The edgeR result for IL1B is therefore NA, not replaced by an unfiltered or modified result.

## Main expression analysis

The complete 63,140-gene matrix was loaded. TMM factors were calculated using all genes, followed by log2-CPM with edgeR's prior count. Differential expression used edgeR quasi-likelihood GLM with `robust=TRUE`, design `~ group`, Control as reference and Sepsis as the coefficient. A Mann–Whitney U test is reported separately and is not mixed with the edgeR P value.

Main results:

- PTGS2: edgeR log2FC -1.0778, P=0.2235, FDR=0.7385; Mann–Whitney P=0.1252.
- IL1B: removed by `filterByExpr`; edgeR log2FC/P/FDR=NA; Mann–Whitney P=0.2769.

## ROC analysis

Sepsis was fixed as the positive class and higher expression was fixed as the predictive direction before calculation. AUC confidence intervals use stratified bootstrap CI with 2,000 resamples (`pROC::ci.auc`). These are apparent discrimination estimates in this 14-sample external cohort, not clinical validation.

- PTGS2: AUC 0.2449, 95% CI 0.0000–0.5510.
- IL1B: AUC 0.6837, 95% CI 0.3469–0.9592.

No joint model was constructed. An AUC below 0.5 for PTGS2 is reported as observed and reflects lower PTGS2 values in Sepsis than Control in this cohort.

Figure 4H was subsequently redrawn from the saved ROC coordinates only. The redraw sorts FPR values and applies a cumulative-monotone TPR envelope before plotting, eliminating threshold-order backtracking and the visually misleading downward vertical segments. AUC values were not recalculated. The redraw script is `work/scripts/redraw_Figure4H_PTGS2_IL1B.R`; identical corrected files are provided as both `Figure_4H_PTGS2_IL1B_ROC.*` and `Figure_4F_PTGS2_IL1B_ROC.*` for panel naming.

## Sensitivity analysis

The primary analysis retains all 14 samples. Secondary analyses separately exclude SRR25368433, separately exclude SRR25368428, and exclude both. These are technical robustness checks only and do not replace the primary analysis. Results are in `differential_sensitivity.tsv` and `ROC_AUC_sensitivity.tsv`.

## Reproducibility

Run:

```bash
Rscript work/scripts/PTGS2_IL1B_edgeR_ROC_analysis.R
```

Software recorded in `sessionInfo.txt` and `package_versions.txt`. The source count matrix, annotation, sample manifest, and earlier alignment/count logs remain in the project tree unchanged.

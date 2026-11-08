# GSE237861 human lung RNA-seq reanalysis

Public reproducibility package for the 14-sample human lung analysis used for the Pharmaceuticals revision (PTGS2 and IL1B). This directory contains the public count-level data, metadata, QC summaries, result tables, figures, and scripts. It deliberately excludes SRA/FASTQ files, BAM/SAM files, HISAT2 indices, reference FASTA/GTF files, caches, and private paths or credentials.

## Reproduce the published count-level results

Requirements: R (4.x), Bioconductor `edgeR`, and CRAN packages `ggplot2` and `pROC`.

From this directory run:

```bash
Rscript scripts/06_downstream/PTGS2_IL1B_edgeR_ROC_analysis.R
Rscript scripts/06_downstream/redraw_Figure4G_PTGS2_IL1B.R
Rscript scripts/06_downstream/redraw_Figure4H_PTGS2_IL1B.R
```

The analysis script reads `data/GSE237861_lung14_gene_counts.tsv`, `data/GSE237861_Ensembl112_gene_annotation.tsv`, and `metadata/lung_14sample_manifest.tsv`, then writes results under `results/`. The original analysis used all genes for `edgeR` filtering/TMM normalization, `-s 1` featureCounts, and Sepsis as ROC-positive. The ROC plotting redraw is plot-only and uses saved coordinates; it removes non-monotone visual connector artifacts without recalculating AUC.

The source project was run on macOS in conda environment `gse237861_rnaseq`. Exact package versions are recorded in `results/package_versions.txt` and `results/sessionInfo.txt`. The public script is the portable counterpart of the executed project script; the source project remains unchanged.

## Data and reference provenance

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237861
- NCBI SRA run accessions and GSM/group assignments: `metadata/lung_14sample_manifest.tsv`
- Ensembl release 112 directory: https://ftp.ensembl.org/pub/release-112/
- GRCh38 primary assembly FASTA: https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
- Ensembl 112 GTF: https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
- Official checksums/file manifest: `docs/ENSEMBL112_primary_assembly_download_manifest.tsv`

The analyzed reference was Ensembl Release 112, GRCh38.p14. The verified source-project SHA-256 values are documented in `docs/reference_download_20261107.log` and the manifest/notes. The large reference files are not redistributed here.

## Analysis order

1. `scripts/01_acquisition/` — acquisition provenance and accessions (raw data are downloaded from NCBI by the researcher).
2. `scripts/02_qc/` — FastQC/MultiQC and fastp provenance; public QC summaries and sanitized fastp JSON reports are in `qc/`.
3. `scripts/03_reference/` — Ensembl 112 reference validation/merge provenance; no reference sequence is bundled.
4. `scripts/04_alignment/` — HISAT2/samtools/featureCounts scripts and original run provenance. These require locally downloaded FASTQ, reference files, and the conda environment.
5. `scripts/05_counting/` — count-matrix construction.
6. `scripts/06_downstream/` — edgeR normalization/differential analysis, ROC analysis, and Figure 4G/H export.

## Results

The full matrix is 63,140 genes × 14 samples. PTGS2 is ENSG00000073756 and IL1B is ENSG00000125538; both mappings are retained in `results/target_gene_annotation.tsv`. Main apparent-discrimination ROC results are PTGS2 AUC 0.245 and IL1B AUC 0.684. PTGS2's AUC below 0.5, IL1B filtering/low expression behavior, non-significant results, and low-depth sensitivity analyses are retained as observed; no sample was removed from the primary analysis.

## Reproducibility status

The downstream count-to-figure workflow is independently rerun from the bundled matrix and metadata and compared with the source results. The upstream FASTQ-to-count workflow is documented with the actual executed scripts/logs, but cannot be end-to-end rerun from this public folder because raw FASTQ, reference files, BAMs, and indices are intentionally excluded. The copied upstream shell scripts retain the historical commands; paths must be supplied through the local project environment before execution.

## License and data use

This package is a reproducibility supplement. Users must follow the GEO/SRA and Ensembl terms of use and cite GSE237861, NCBI SRA, Ensembl, HISAT2, samtools, featureCounts/Subread, edgeR, and pROC as appropriate.

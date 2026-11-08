#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
BIN="<CONDA_ENV_BIN>"
FA="$PROJECT/work/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="$PROJECT/work/reference/browser_downloads_20261107/Homo_sapiens.GRCh38.112.gtf"
IDX="$PROJECT/work/reference/hisat2_index/GRCh38_Ensembl112"
FQ="$PROJECT/work/trimmed/SRR25368441.trimmed.fastq.gz"

# Index construction (already completed for this project):
# "$BIN/hisat2-build" -p 8 "$FA" "$IDX"

# Single-end QuantSeq FWD test alignment. FWD reads directly reflect mRNA sequence;
# --rna-strandness F is recorded for the HISAT2 alignment run.
mkdir -p "$PROJECT/work/aligned" "$PROJECT/work/counts" "$PROJECT/work/logs"
"$BIN/hisat2" -p 8 --dta --rna-strandness F -x "$IDX" -U "$FQ" \
  --summary-file "$PROJECT/work/logs/SRR25368441.hisat2.summary.txt" \
  2> "$PROJECT/work/logs/SRR25368441.hisat2_alignment.stderr.log" \
  | "$BIN/samtools" sort -@ 4 -m 1G -o "$PROJECT/work/aligned/SRR25368441.hisat2.sorted.bam" -
"$BIN/samtools" index -@ 4 "$PROJECT/work/aligned/SRR25368441.hisat2.sorted.bam"
"$BIN/samtools" quickcheck -v "$PROJECT/work/aligned/SRR25368441.hisat2.sorted.bam"
"$BIN/samtools" flagstat "$PROJECT/work/aligned/SRR25368441.hisat2.sorted.bam" > "$PROJECT/work/logs/SRR25368441.flagstat.txt"

# Single-end, exon-level, gene_id aggregation; test all strand settings.
for s in 0 1 2; do
  "$BIN/featureCounts" -T 4 -a "$GTF" -o "$PROJECT/work/counts/SRR25368441.s${s}.counts.txt" \
    -t exon -g gene_id -s "$s" "$PROJECT/work/aligned/SRR25368441.hisat2.sorted.bam" \
    > "$PROJECT/work/logs/SRR25368441_featureCounts_s${s}.log" 2>&1
done

#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
BIN="${BIN_DIR:-$(dirname "$(command -v hisat2)")}"
IDX="$PROJECT/work/reference/hisat2_index/GRCh38_Ensembl112"
GTF="$PROJECT/work/reference/browser_downloads_20261107/Homo_sapiens.GRCh38.112.gtf"
MANIFEST="$PROJECT/work/metadata/lung_14sample_manifest.tsv"

mkdir -p "$PROJECT/work/aligned" "$PROJECT/work/counts" "$PROJECT/work/logs"
QC="$PROJECT/work/logs/all14_alignment_count_qc.tsv"
if [[ ! -s "$QC" ]]; then
  printf 'GSM\tSRX\tSRR\tgroup\ttrimmed_reads\thisat2_unique_reads\thisat2_multimapped_reads\thisat2_unmapped_reads\thisat2_overall_pct\tprimary_mapped_records\tsecondary_records\tfeatureCounts_total_records\tfeatureCounts_assigned\tfeatureCounts_assigned_pct\tfeatureCounts_nofeatures\tfeatureCounts_ambiguity\tBAM_bytes\n' > "$QC"
fi

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r gsm srx srr group tissue title; do
  fq="$PROJECT/work/trimmed/${srr}.trimmed.fastq.gz"
  bam="$PROJECT/work/aligned/${srr}.hisat2.sorted.bam"
  summary="$PROJECT/work/logs/${srr}.hisat2.summary.txt"
  stderr="$PROJECT/work/logs/${srr}.hisat2_alignment.stderr.log"
  flagstat="$PROJECT/work/logs/${srr}.flagstat.txt"
  count="$PROJECT/work/counts/${srr}.s1.counts.txt"
  countlog="$PROJECT/work/logs/${srr}_featureCounts_s1.log"
  test -s "$fq"

  if [[ "$srr" != "SRR25368441" && ! -s "$bam" ]]; then
    "$BIN/hisat2" -p 8 --dta --rna-strandness F -x "$IDX" -U "$fq" \
      --summary-file "$summary" 2> "$stderr" \
      | "$BIN/samtools" sort -@ 4 -m 1G -o "$bam" -
    "$BIN/samtools" index -@ 4 "$bam"
    "$BIN/samtools" quickcheck -v "$bam"
    "$BIN/samtools" flagstat "$bam" > "$flagstat"
  fi

  if [[ "$srr" != "SRR25368441" && ! -s "$count" ]]; then
    "$BIN/featureCounts" -T 4 -a "$GTF" -o "$count" -t exon -g gene_id -s 1 "$bam" > "$countlog" 2>&1
  fi

  reads=$(gzip -cd "$fq" | awk 'END{print NR/4}')
  unique=$(awk '/aligned exactly 1 time/{print $1; exit}' "$summary")
  multi=$(awk '/aligned >1 times/{print $1; exit}' "$summary")
  unmapped=$(awk '/aligned 0 times/{print $1; exit}' "$summary")
  overall=$(awk '/overall alignment rate/{print $1; exit}' "$summary")
  primary_mapped=$(awk '/primary mapped/{print $1; exit}' "$flagstat")
  secondary=$(awk '/secondary/{print $1; exit}' "$flagstat")
  total=$(awk -F '\t' '$1=="Assigned"{a=$2} $1=="Unassigned_Unmapped"{b=$2} $1=="Unassigned_MultiMapping"{c=$2} $1=="Unassigned_NoFeatures"{d=$2} $1=="Unassigned_Ambiguity"{e=$2} END{print a+b+c+d+e}' "${count}.summary")
  assigned=$(awk -F '\t' '$1=="Assigned"{print $2}' "${count}.summary")
  nofeatures=$(awk -F '\t' '$1=="Unassigned_NoFeatures"{print $2}' "${count}.summary")
  ambiguity=$(awk -F '\t' '$1=="Unassigned_Ambiguity"{print $2}' "${count}.summary")
  assigned_pct=$(awk -v a="$assigned" -v t="$total" 'BEGIN{printf "%.2f",100*a/t}')
  bam_bytes=$(stat -f %z "$bam")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$gsm" "$srx" "$srr" "$group" "$reads" "$unique" "$multi" "$unmapped" "$overall" "$primary_mapped" "$secondary" "$total" "$assigned" "$assigned_pct" "$nofeatures" "$ambiguity" "$bam_bytes" >> "$QC"
done

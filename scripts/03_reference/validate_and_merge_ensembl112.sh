#!/usr/bin/env bash
set -euo pipefail

# Reproducible reference preparation for Ensembl Release 112 / GRCh38.p14.
# This script intentionally does not download, delete, or build HISAT2 indexes.
PROJECT="${PROJECT_ROOT:-$(pwd)}"
REF="$PROJECT/work/reference"
DOWNLOADS="$HOME/Downloads"
BROWSER="$REF/browser_downloads_20261107"
OUT="$REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

mkdir -p "$REF/validation_20261107"
: > "$OUT"
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
  if [[ "$n" == 1 ]]; then f="$BROWSER/Homo_sapiens.GRCh38.dna.chromosome.1.fa"; else f="$DOWNLOADS/Homo_sapiens.GRCh38.dna.chromosome.$n.fa"; fi
  test -s "$f"
  cat "$f" >> "$OUT"
done
cat "$DOWNLOADS/Homo_sapiens.GRCh38.dna.nonchromosomal.fa" >> "$OUT"

awk 'BEGIN{h=0;len=0;bad=0;dup=0;total=0} /^>/{if(h>0){if(len==0)bad++;total+=len};h++;id=$1;sub(/^>/,"",id);if(seen[id]++)dup++;len=0;next}{gsub(/[[:space:]]/,"");if($0~/[^ACGTNacgtn]/)bad++;len+=length}END{total+=len;print "merged_sequences",h,"duplicate_ids",dup,"invalid_or_empty",bad,"total_bases",total}' "$OUT" > "$REF/validation_20261107/merged_fasta_summary.txt"

GTF="$BROWSER/Homo_sapiens.GRCh38.112.gtf"
awk 'BEGIN{FS="\t";bad=0;data=0;comments=0}/^#/{comments++;next}{data++;if(NF!=9)bad++;if($4!~/^[0-9]+$/||$5!~/^[0-9]+$/||$4>$5)bad++}END{print "gtf_data_lines",data,"comments",comments,"malformed",bad}' "$GTF" > "$REF/validation_20261107/gtf_structure_summary.txt"
shasum -a 256 "$OUT" > "$REF/validation_20261107/primary_assembly.fa.sha256"
shasum -a 256 "$GTF" > "$REF/validation_20261107/Homo_sapiens.GRCh38.112.gtf.sha256"

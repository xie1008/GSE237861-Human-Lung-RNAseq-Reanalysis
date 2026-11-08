#!/usr/bin/env python3
import csv, re, os
from pathlib import Path

PROJECT = Path(os.environ.get('PROJECT_ROOT', Path.cwd()))
manifest = PROJECT/'work/metadata/lung_14sample_manifest.tsv'
gtf = PROJECT/'work/reference/browser_downloads_20261107/Homo_sapiens.GRCh38.112.gtf'
count_dir = PROJECT/'work/counts'
out_dir = count_dir/'matrix'
out_dir.mkdir(parents=True, exist_ok=True)

rows = list(csv.DictReader(manifest.open(), delimiter='\t'))
samples = [r['SRR'] for r in rows]
gsm = {r['SRR']: r['GSM'] for r in rows}

counts = {}
gene_order = None
for srr in samples:
    path = count_dir/f'{srr}.s1.counts.txt'
    d = {}
    with path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] == 'Geneid':
                continue
            d[fields[0]] = int(fields[-1])
    if gene_order is None:
        gene_order = list(d)
    elif list(d) != gene_order:
        raise RuntimeError(f'gene ID/order mismatch in {srr}')
    counts[srr] = d

matrix = out_dir/'GSE237861_lung14_gene_counts.tsv'
with matrix.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['gene_id'] + [gsm[s] for s in samples])
    for gene in gene_order:
        w.writerow([gene] + [counts[s][gene] for s in samples])

attrs_re = re.compile(r'(gene_id|gene_name|gene_biotype) "([^"]+)"')
ann = {}
with gtf.open() as f:
    for line in f:
        if line.startswith('#'):
            continue
        x = line.rstrip('\n').split('\t')
        if len(x) != 9 or x[2] != 'gene':
            continue
        a = dict(attrs_re.findall(x[8]))
        gid = a.get('gene_id')
        if gid and gid not in ann:
            ann[gid] = (a.get('gene_name',''), a.get('gene_biotype',''), x[0], x[6])

annotation = out_dir/'GSE237861_Ensembl112_gene_annotation.tsv'
with annotation.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['gene_id','gene_name','gene_biotype','seqname','strand'])
    for gene in gene_order:
        w.writerow([gene] + list(ann.get(gene, ('','','',''))))

targets = {'PTGS2','PTGES','STAT3','ALOX5'}
target_path = out_dir/'target_gene_mapping.tsv'
with target_path.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['gene_id','gene_name','gene_biotype','seqname','strand'])
    for gene in gene_order:
        if ann.get(gene, ('','','',''))[0] in targets:
            w.writerow([gene] + list(ann[gene]))

print('samples', len(samples), 'genes', len(gene_order), 'matrix', matrix)
print('target rows', sum(1 for gene in gene_order if ann.get(gene, ('','','',''))[0] in targets))

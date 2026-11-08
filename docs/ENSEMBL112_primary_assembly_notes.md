# Ensembl Release 112 GRCh38 primary assembly plan

Checked against the official Ensembl DNA directory and its `README`:

- Assembly accession: `GCA_000001405.29` (GRCh38.p14).
- Ensembl defines primary assembly as all top-level sequence regions excluding haplotypes and patches.
- Therefore primary assembly is not just the 25 chromosome files. It includes the 25 chromosome sequences plus the nonchromosomal/unassigned top-level sequences in `Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz`.
- The official assembly inventory reports 194 primary top-level sequences: 25 chromosome sequences and 169 scaffold/nonchromosomal sequences, totaling 3,099,750,718 bp. The nonchromosomal component contributes 11,464,317 bp. These counts and lengths were obtained from the Ensembl assembly endpoint for the same assembly accession and are to be rechecked against the downloaded FASTA headers before final promotion.
- `dna.alt.fa.gz`, `dna.patch.fa.gz`, haplotype files, and `dna.toplevel.fa.gz` must not be added to the primary-assembly FASTA. `toplevel` includes haplotype/patch regions; primary assembly explicitly excludes them.
- The uncompressed browser files already present for chromosome 1, chromosome 2, and the Release 112 GTF are retained and must not be downloaded again.

Official directory:

https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/

Official README:

https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/README

Full primary assembly (do not download again unless explicitly needed):

https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Required browser downloads are listed in `ENSEMBL112_primary_assembly_download_manifest.tsv`. The missing component is the nonchromosomal file plus chromosome 3–22, X, Y, and MT; chr1, chr2, and the GTF are already present.

Space plan: the uncompressed primary assembly is approximately 3.10 GB based on the official sequence inventory. Keep one download set and one merged FASTA only; reserve additional space for the future HISAT2 index and temporary index files. With approximately 41 GiB currently free, the reference preparation is feasible, but index construction should wait until the download set is complete and temporary files are cleaned or moved deliberately.

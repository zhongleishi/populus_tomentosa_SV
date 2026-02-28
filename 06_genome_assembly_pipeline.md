# 06. Genome Assembly Pipeline (PacBio + Illumina + 10X + Hi-C)

This document provides  a complete runnable bash script (`run_assembly.sh`) embedded below.

The script performs:
1) k-mer survey (Jellyfish + GenomeScope input)
2) primary assembly (Canu)
3) contiguity improvement (HERA) with default parameters
4) redundancy removal (Redundans)
5) haplotig purging (Purge Haplotigs)
6) deterministic self-alignment-based merging (minimap2 all-vs-all + keep-longest-per-cluster)

---

## Software
- Jellyfish
- GenomeScope (manual step using histogram)
- Canu
- HERA
- Redundans
- Purge Haplotigs
- minimap2, samtools
- seqtk (for contig extraction)

---

## Inputs (edit in script)
- Illumina: `data/illumina_R1.fastq.gz`, `data/illumina_R2.fastq.gz`
- PacBio: `data/pacbio.fastq.gz`
- 10X: `data/10x_R1.fastq.gz`, `data/10x_R2.fastq.gz`

---

## Outputs
- Final assembly: `06_merge/final_assembly.fasta`
- Intermediate directories: `01_kmer_analysis/` ... `06_merge/`

---

## How to run

1) Copy the script block below into a file named `run_assembly.sh`
2) Run:

```bash
bash run_assembly.sh 2>&1 | tee assembly.log
```

---

## Script: run_assembly.sh

```bash
#!/bin/bash
set -euo pipefail

THREADS=40
GENOME_SIZE="390m"

ILLUMINA_R1="data/illumina_R1.fastq.gz"
ILLUMINA_R2="data/illumina_R2.fastq.gz"
PACBIO_READS="data/pacbio.fastq.gz"
TENX_R1="data/10x_R1.fastq.gz"
TENX_R2="data/10x_R2.fastq.gz"

OUT_KMER="01_kmer_analysis"
OUT_CANU="02_canu_assembly"
OUT_HERA="03_hera_assembly"
OUT_REDUNDANS="04_redundans"
OUT_PURGE="05_purge_haplotigs"
OUT_MERGE="06_merge"

mkdir -p "$OUT_KMER" "$OUT_CANU" "$OUT_HERA" "$OUT_REDUNDANS" "$OUT_PURGE" "$OUT_MERGE"

echo "[1/6] K-mer survey (Jellyfish + GenomeScope input)"
jellyfish count -C -m 21 -s 10G -t "$THREADS" -o "$OUT_KMER/reads.jf" <(zcat "$ILLUMINA_R1" "$ILLUMINA_R2")
jellyfish histo -t "$THREADS" "$OUT_KMER/reads.jf" > "$OUT_KMER/reads.histo"
echo "GenomeScope input: $OUT_KMER/reads.histo"

echo "[2/6] Primary assembly (Canu)"
canu -p 1316 -d "$OUT_CANU" \
  genomeSize="$GENOME_SIZE" \
  useGrid=false \
  maxThreads="$THREADS" \
  -pacbio "$PACBIO_READS"

DRAFT_CONTIGS="$OUT_CANU/1316.contigs.fasta"
test -s "$DRAFT_CONTIGS"

echo "[3/6] Contiguity improvement (HERA, default parameters)"
# Default parameters used for reproducibility:
# - min overlap length: 5000 bp
# - min identity: 0.97
# If your HERA CLI differs, adapt the command to your installation while keeping equivalent thresholds.
cp "$DRAFT_CONTIGS" "$OUT_HERA/draft_assembly.fasta"

HERA_INPUT="$OUT_HERA/draft_assembly.fasta"
HERA_OUT="$OUT_HERA/hera.fasta"

HERA.pl \
  --genome "$HERA_INPUT" \
  --reads "$PACBIO_READS" \
  --min_ovlp 5000 \
  --min_id 0.97 \
  --threads "$THREADS" \
  --out "$HERA_OUT"

test -s "$HERA_OUT"

echo "[4/6] Redundancy removal (Redundans)"
redundans.py -v \
  -i "$ILLUMINA_R1" "$ILLUMINA_R2" \
  -f "$HERA_OUT" \
  -o "$OUT_REDUNDANS" \
  -t 10 \
  -identity 0.55 \
  -overlap 0.80 \
  --noscaffolding \
  --nogapclosing

CURRENT_DRAFT="$OUT_REDUNDANS/scaffolds.filled.fa"
test -s "$CURRENT_DRAFT"

echo "[5/6] Haplotig purging (Purge Haplotigs)"
minimap2 -t "$THREADS" -ax map-pb "$CURRENT_DRAFT" "$PACBIO_READS" \
  | samtools sort -@ "$THREADS" -o "$OUT_PURGE/aligned.bam" -
samtools index "$OUT_PURGE/aligned.bam"

purge_haplotigs readhist -b "$OUT_PURGE/aligned.bam" -g "$CURRENT_DRAFT" -t "$THREADS"

# Default thresholds (from the original notes; adjust if histogram suggests different cutoffs)
purge_haplotigs contigcov -i "$OUT_PURGE/aligned.bam.gencov" -l 5 -m 30 -h 150
purge_haplotigs purge -g "$CURRENT_DRAFT" -c "$OUT_PURGE/coverage_stats.csv" -t "$THREADS" -o "$OUT_PURGE/curated"

PURGED="$OUT_PURGE/curated.fasta"
test -s "$PURGED"

echo "[6/6] Deterministic self-alignment merging (minimap2 asm5 + longest-per-cluster)"
minimap2 -x asm5 -t "$THREADS" "$PURGED" "$PURGED" > "$OUT_MERGE/overlaps.paf"

python - << 'PY'
import sys
from collections import defaultdict

paf = "06_merge/overlaps.paf"
fa = "05_purge_haplotigs/curated.fasta"
out = "06_merge/keep_ids.txt"

# fasta lengths
lengths = {}
with open(fa) as f:
    name = None
    L = 0
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                lengths[name] = L
            name = line[1:].split()[0]
            L = 0
        else:
            L += len(line)
    if name is not None:
        lengths[name] = L

parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

# thresholds (paper-default)
MIN_ID = 0.97
MIN_COV = 0.80

with open(paf) as f:
    for line in f:
        if not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        q = p[0]; qlen = int(p[1])
        t = p[5]
        if q == t:
            continue
        nmatch = int(p[9]); alnlen = int(p[10])
        if alnlen == 0 or qlen == 0:
            continue
        identity = nmatch / alnlen
        cov = alnlen / qlen
        if identity >= MIN_ID and cov >= MIN_COV:
            union(q, t)

clusters = defaultdict(list)
for cid in lengths:
    clusters[find(cid)].append(cid)

keep = []
for root, items in clusters.items():
    items.sort(key=lambda x: lengths.get(x, 0), reverse=True)
    keep.append(items[0])

with open(out, "w") as w:
    for k in sorted(keep):
        w.write(k + "\n")

print(f"Clusters: {len(clusters)}; kept contigs: {len(keep)}", file=sys.stderr)
PY

seqtk subseq "$PURGED" "$OUT_MERGE/keep_ids.txt" > "$OUT_MERGE/final_assembly.fasta"
echo "Final assembly: $OUT_MERGE/final_assembly.fasta"
```

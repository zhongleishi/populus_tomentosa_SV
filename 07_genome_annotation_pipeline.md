# 07. Genome Annotation Pipeline (Trinity + HISAT2/StringTie + PASA + TransDecoder + Augustus + EVM)

This document provides a  complete runnable bash script (`run_annotation.sh`) embedded below.

The script performs:
1) transcriptome assembly: **Trinity** (de novo) + **HISAT2/StringTie** (reference-guided)
2) transcript alignment & gene structure modeling: **PASA**
3) ORF prediction: **TransDecoder**
4) build a deterministic PASA-derived training gene set for Augustus (default rules)
5) train and predict genes: **Augustus (autoAug)**
6) integrate evidence into consensus gene models: **EVM**

---

## Software (recommended versions)
- Trinity (>=2.8)
- HISAT2, StringTie
- PASA (>=2.5)
- TransDecoder
- Augustus (+ autoAug.pl)
- EvidenceModeler (EVM)
- samtools, gffread
- MySQL (for PASA)

---

## Inputs (edit in script)
- Genome FASTA: `data/genome_1316.fasta`
- RNA-seq reads: `data/rna_1.fastq.gz`, `data/rna_2.fastq.gz`
- PASA config: `conf/alignAssembly.config` (must define DB connection and DB name)

---

## Outputs
- Consensus annotation: `06_evm/evm_consensus.gff3`

---

## How to run

1) Copy the script block below into a file named `run_annotation.sh`
2) Run:

```bash
bash run_annotation.sh 2>&1 | tee annotation.log
```

---

## Script: run_annotation.sh

```bash
#!/bin/bash
set -euo pipefail

THREADS=40
GENOME_FASTA="data/genome_1316.fasta"
RNA_SEQ_LEFT="data/rna_1.fastq.gz"
RNA_SEQ_RIGHT="data/rna_2.fastq.gz"

OUT_TRINITY="01_trinity_denovo"
OUT_REF="02_ref_guided"
OUT_PASA="03_pasa"
OUT_TD="04_transdecoder"
OUT_AUG="05_augustus"
OUT_EVM="06_evm"

mkdir -p "$OUT_TRINITY" "$OUT_REF" "$OUT_PASA" "$OUT_TD" "$OUT_AUG" "$OUT_EVM"

echo "[1/6] Trinity"
Trinity --seqType fq --max_memory 50G --CPU "$THREADS" \
  --left "$RNA_SEQ_LEFT" --right "$RNA_SEQ_RIGHT" \
  --output "$OUT_TRINITY"

echo "[2/6] HISAT2 + StringTie"
hisat2-build "$GENOME_FASTA" "$OUT_REF/genome_index"

hisat2 -p "$THREADS" -x "$OUT_REF/genome_index" \
  -1 "$RNA_SEQ_LEFT" -2 "$RNA_SEQ_RIGHT" \
  -S "$OUT_REF/aligned.sam"

samtools sort -@ "$THREADS" -o "$OUT_REF/aligned.bam" "$OUT_REF/aligned.sam"
samtools index "$OUT_REF/aligned.bam"

stringtie "$OUT_REF/aligned.bam" -p "$THREADS" -o "$OUT_REF/transcripts.gtf"
gffread "$OUT_REF/transcripts.gtf" -g "$GENOME_FASTA" -w "$OUT_REF/transcripts.fasta"

echo "[3/6] PASA"
cat "$OUT_TRINITY/Trinity.fasta" "$OUT_REF/transcripts.fasta" > "$OUT_PASA/all_transcripts.fasta"
seqclean "$OUT_PASA/all_transcripts.fasta"

cp conf/alignAssembly.config "$OUT_PASA/"

Launch_PASA_pipeline.pl \
  -c "$OUT_PASA/alignAssembly.config" \
  -C -R -g "$GENOME_FASTA" \
  -t "$OUT_PASA/all_transcripts.fasta.clean" \
  --ALIGNERS blat,gmap --CPU "$THREADS"

echo "[4/6] TransDecoder"
# Auto-detect PASA output names (DB-dependent)
PASA_ASM_FASTA=$(ls "$OUT_PASA"/*.assemblies.fasta | head -n 1)
PASA_ASM_GFF3=$(ls "$OUT_PASA"/*.pasa_assemblies.gff3 | head -n 1)

test -s "$PASA_ASM_FASTA"
test -s "$PASA_ASM_GFF3"

cp "$PASA_ASM_FASTA" "$OUT_TD/pasa_assemblies.fasta"

TransDecoder.LongOrfs -t "$OUT_TD/pasa_assemblies.fasta"
TransDecoder.Predict -t "$OUT_TD/pasa_assemblies.fasta"

echo "[5/6] Build PASA training gene set (default deterministic rules)"
# Default rule for reproducibility:
# - keep transcript models that are multi-exon (>=2 exons) AND contain CDS features in PASA GFF3.
python - << 'PY'
import re
from collections import defaultdict

gff = open("03_pasa/" + __import__("glob").glob("03_pasa/*.pasa_assemblies.gff3")[0].split("/")[-1])
out = open("05_augustus/pasa_training_genes.gff3", "w")

exon_count = defaultdict(int)
has_cds = set()
lines_by_tid = defaultdict(list)

for line in gff:
    if line.startswith("#") or not line.strip():
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9:
        continue
    ftype = p[2].lower()
    attrs = p[8]

    m = re.search(r"(?:Parent|ID)=([^;]+)", attrs)
    if not m:
        continue
    tid = m.group(1)

    lines_by_tid[tid].append(line)
    if ftype == "exon":
        exon_count[tid] += 1
    if ftype == "cds":
        has_cds.add(tid)

kept = [tid for tid, c in exon_count.items() if c >= 2 and tid in has_cds]

for tid in kept:
    for l in lines_by_tid[tid]:
        out.write(l)

print(f"Training transcripts kept: {len(kept)}")
PY

echo "[5/6] Augustus training + prediction"
gff2gbSmallDNA.pl "$OUT_AUG/pasa_training_genes.gff3" "$GENOME_FASTA" 1000 "$OUT_AUG/genes.gb"

autoAug.pl \
  --species=populus_1316 \
  --genome="$GENOME_FASTA" \
  --trainingset="$OUT_AUG/genes.gb" \
  --cdna="$OUT_PASA/all_transcripts.fasta" \
  --singleCPU

augustus --species=populus_1316 "$GENOME_FASTA" > "$OUT_AUG/augustus_preds.gff"

echo "[6/6] EVM"
cat << EOF > "$OUT_EVM/weights.txt"
TRANSCRIPT  PASA  10
ABINITIO_PREDICTION  Augustus  5
EOF

EVidenceModeler \
  --genome "$GENOME_FASTA" \
  --gene_predictions "$OUT_AUG/augustus_preds.gff" \
  --transcript_alignments "$PASA_ASM_GFF3" \
  --weights "$OUT_EVM/weights.txt" \
  --cpu "$THREADS" \
  > "$OUT_EVM/evm_consensus.gff3"

echo "Annotation finished: $OUT_EVM/evm_consensus.gff3"
```

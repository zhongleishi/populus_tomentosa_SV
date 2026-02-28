# 03. RNA-seq Alignment (HISAT2) and BAM generation

## Software
- HISAT2
- SAMtools

## Inputs
- Reference FASTA: `reference.fa`
- Paired reads: `sample1.good.1.fq.gz`, `sample1.good.2.fq.gz`

## Outputs
- HISAT2 index prefix: `ref`
- Sorted BAM: `sorted/sample1.sorted.bam`

---

## Step 1 — Build HISAT2 index

```bash
hisat2-build -f reference.fa -p 7 ref
```

---

## Step 2 — Align

```bash
hisat2 -q -x ref \
  --min-intronlen 20 \
  --max-intronlen 4000 \
  --rna-strandness RF \
  -1 sample1.good.1.fq.gz \
  -2 sample1.good.2.fq.gz \
  -S sample1.sam \
  --dta-cufflinks \
  2> sample1.hisat2.log
```

---

## Step 3 — Sort + index BAM

```bash
mkdir -p sorted
samtools view -Su sample1.sam | samtools sort -@ 4 -o sorted/sample1.sorted.bam
samtools index sorted/sample1.sorted.bam
```

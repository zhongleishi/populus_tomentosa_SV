# 01. Structural Variant (SV) Calling and Population Genotyping

This workflow provides a reproducible SV discovery and population genotyping procedure:
1) per-sample SV calling with **Manta**  
2) per-sample filtering to retain high-confidence (precise) SVs  
3) multi-sample SV merging using **svimmer**  
4) population SV genotyping using **Graphtyper**  
5) merge per-chromosome results and final filtering

---

## Software
- Manta
- svimmer
- Graphtyper
- bcftools
- bgzip/tabix
- awk

---

## Inputs
- Reference FASTA: `genome.fasta`
- BAMs: e.g. `pto1.bam`, `pto2.bam`, ...
- Chromosome IDs: `Chr01..Chr19 ChrX` (adjust to your reference)
- `bam.list`: one BAM path per line (for Graphtyper)

---

## Step 1 — Per-sample SV calling (Manta)

Configure (example with multiple BAMs):

```bash
python configManta.py \
  --referenceFasta genome.fasta \
  --runDir manta_run \
  --bam pto1.bam --bam pto2.bam --bam pto3.bam
```

Run workflow:

```bash
manta_run/runWorkflow.py -m local -j 40
```

Typical output:
- `manta_run/results/variants/diploidSV.vcf.gz`

---

## Step 2 — Quality control / precise SV filtering

Filtering rules:
- keep header (`^#`)
- `FILTER == PASS`
- `QUAL > 20`
- remove records containing `IMPRECISE`

```bash
zcat manta_run/results/variants/diploidSV.vcf.gz \
  | awk '/^#/ || ($7 == "PASS" && $6 > 20 && $8 !~ /IMPRECISE/)' \
  | bgzip -c > pto1.precise.vcf.gz

tabix -p vcf pto1.precise.vcf.gz
```

Repeat for each sample VCF.

---

## Step 3 — Merge SVs across individuals (svimmer)

Create `vcf.list` (one VCF per line), e.g.:
```text
pto1.precise.vcf.gz
pto2.precise.vcf.gz
pto3.precise.vcf.gz
```

Run (example chromosome list; adjust as needed):

```bash
svimmer --threads 60 vcf.list Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 ChrX
```

---

## Step 4 — Population genotyping (Graphtyper)

Genotype SVs per chromosome (example Chr01):

```bash
graphtyper genotype_sv genome.fasta sample_merged_sort.vcf.gz \
  --sams=bam.list \
  --region Chr01 \
  --threads 40
```

---

## Step 5 — Merge per-chromosome Graphtyper results + final filter

Collect all VCFs:

```bash
echo Chr{01..19} ChrX | tr ' ' '\n' | while read chrom; do
  if [[ ! -d sv_results/${chrom} ]]; then
    continue
  fi
  find sv_results/${chrom} -name "*.vcf.gz" | sort
done > vcf_file_list
```

Concatenate:

```bash
bcftools concat --naive --file-list vcf_file_list -Oz -o graphtyper.all.vcf.gz
tabix -p vcf graphtyper.all.vcf.gz
```

PASS-only:

```bash
bcftools view -f PASS -Oz -o graphtyper.PASS.vcf.gz graphtyper.all.vcf.gz
tabix -p vcf graphtyper.PASS.vcf.gz
```

Optional SV-type filter (adjust to your VCF encoding):

```bash
zcat graphtyper.PASS.vcf.gz \
  | awk '/^#/ || ($5 ~ /AGGREGATED|INV/)' \
  | bgzip -c > graphtyper.PASS.filtered.vcf.gz
tabix -p vcf graphtyper.PASS.filtered.vcf.gz
```

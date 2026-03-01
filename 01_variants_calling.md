# Variants Calling (SV + SNP/INDEL) Workflow

This document provides a reproducible variant discovery workflow including:

- **SV calling & population genotyping**
  1) per-sample SV calling with **Manta**  
  2) per-sample filtering to retain high-confidence (precise) SVs  
  3) multi-sample SV merging using **svimmer**  
  4) population SV genotyping using **Graphtyper**  
  5) merge per-chromosome results and final filtering

- **SNP/INDEL calling (GATK joint genotyping)**
  1) reference preparation  
  2) read mapping (**BWA-MEM**) + BAM sort  
  3) mark duplicates (**GATK MarkDuplicates**) + index  
  4) per-sample variant calling (**GATK HaplotypeCaller** in **GVCF** mode)  
  5) merge GVCFs (**CombineGVCFs**)  
  6) joint genotyping (**GenotypeGVCFs**)  
  7) split SNP vs INDEL (**SelectVariants**)  
  8) hard filtering (**VariantFiltration**)

---

## Software

### SV
- Manta
- svimmer
- Graphtyper
- bcftools
- bgzip/tabix
- awk

### SNP/INDEL
- bwa
- samtools
- GATK (4.x)
- bgzip/tabix (optional, for indexing vcf.gz)

---

## Inputs

### Common
- Reference FASTA: `genome.fasta` (example name; adjust)
- Chromosome IDs: `Chr01..Chr19 ChrX` (adjust to your reference)

### SV
- BAMs: e.g. `pto1.bam`, `pto2.bam`, ...
- `bam.list`: one BAM path per line (for Graphtyper)

### SNP/INDEL
- Clean paired FASTQs per sample: `sample.1.clean.fq.gz` and `sample.2.clean.fq.gz`

---

# Part A. Structural Variant (SV) Calling and Population Genotyping

## Step A1 — Per-sample SV calling (Manta)

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

## Step A2 — Quality control / precise SV filtering

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

## Step A3 — Merge SVs across individuals (svimmer)

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

## Step A4 — Population genotyping (Graphtyper)

Genotype SVs per chromosome (example Chr01):

```bash
graphtyper genotype_sv genome.fasta sample_merged_sort.vcf.gz \
  --sams=bam.list \
  --region Chr01 \
  --threads 40
```

---

## Step A5 — Merge per-chromosome Graphtyper results + final filter

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

---

# Part B. SNP/INDEL Calling (GATK joint genotyping)

## Step B0 — Reference preparation (建库)

> Note: for GATK, the dictionary file name is typically derived from the FASTA prefix.
> If your reference is `07.1316.genome.sorted.fasta`, the dict should be `07.1316.genome.sorted.dict`.

```bash
# BWA index
bwa index ref.fa

# samtools faidx
samtools faidx ref.fa

# GATK sequence dictionary
gatk CreateSequenceDictionary -R ref.fa -O ref.dict
```

---

## Step B1 — Mapping (BWA-MEM) + sort BAM

Input naming convention:
- `sample.1.clean.fq.gz`
- `sample.2.clean.fq.gz`

Generate mapping commands:

```bash
for i in `ls *.1.clean.fq.gz`; do
  m=${i/.1.clean.fq.gz/}

  # BWA read group: ID/SM use sample name; PL fixed as illumina
  echo "bwa mem -t 10 -R \"@RG\\tID:${m}\\tPL:illumina\\tSM:${m}\" ref.fa ${m}.1.clean.fq.gz ${m}.2.clean.fq.gz \
    2> ${m}.bwa.log \
    | samtools view -bS -q 20 \
    | samtools sort -@ 5 -o ${m}.sorted.bam -"

  echo "samtools flagstat ${m}.sorted.bam > ${m}.sorted.bam.flagstat.txt"
done > mapping1.list.sh
```

Run:
```bash
bash mapping1.list.sh
```

---

## Step B2 — MarkDuplicates + index

```bash
for i in `ls *.sorted.bam`; do
  m=${i/.sorted.bam/}
  echo "gatk MarkDuplicates -I ${m}.sorted.bam -O ${m}.sorted.markdup.bam -M ${m}.sorted.markdup_metrics.txt"
  echo "samtools index ${m}.sorted.markdup.bam"
done > gatk.md.list.sh
```

Run:
```bash
bash gatk.md.list.sh
```

---

## Step B3 — Per-sample variant calling (HaplotypeCaller in GVCF mode)

```bash
for i in `ls *.sorted.markdup.bam`; do
  m=${i/.sorted.markdup.bam/}

  echo "gatk --java-options \"-Xmx40g -XX:ParallelGCThreads=4\" HaplotypeCaller \
    -R ref.fa \
    -I ${i} \
    --emit-ref-confidence GVCF \
    -O ${m}.g.vcf.gz"
done > gvcf.list.sh
```

Run:
```bash
bash gvcf.list.sh
```

---

## Step B4 — Combine GVCFs (merge g.vcf files)

Create a text file listing your per-sample GVCFs (one per line), e.g. `gvcf.list`:
```text
pto1.g.vcf.gz
pto2.g.vcf.gz
pto3.g.vcf.gz
```

Then run:

```bash
gatk CombineGVCFs \
  -R ref.fa \
  --variant gvcf.list \
  -O cohort.g.vcf.gz
```

---

## Step B5 — Joint genotyping (GenotypeGVCFs)

```bash
gatk GenotypeGVCFs \
  -R ref.fa \
  -V cohort.g.vcf.gz \
  -O cohort.vcf.gz
```

---

## Step B6 — Split SNP and INDEL (SelectVariants)

```bash
# SNPs
gatk SelectVariants \
  -V cohort.vcf.gz \
  -O cohort.snp.vcf.gz \
  --select-type-to-include SNP

# INDELs
gatk SelectVariants \
  -V cohort.vcf.gz \
  -O cohort.indel.vcf.gz \
  --select-type-to-include INDEL
```

---

## Step B7 — Hard filtering (VariantFiltration)

### SNP hard filter

```bash
gatk VariantFiltration \
  -V cohort.snp.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O cohort.snp.hardfilter.vcf.gz
```

### INDEL hard filter

```bash
gatk VariantFiltration \
  -V cohort.indel.vcf.gz \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "SOR > 10.0" --filter-name "SOR10" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O cohort.indel.hardfilter.vcf.gz
```

### (Optional) extract PASS-only variants

```bash
bcftools view -f PASS -Oz -o cohort.snp.PASS.vcf.gz cohort.snp.hardfilter.vcf.gz
bcftools view -f PASS -Oz -o cohort.indel.PASS.vcf.gz cohort.indel.hardfilter.vcf.gz

tabix -p vcf cohort.snp.PASS.vcf.gz
tabix -p vcf cohort.indel.PASS.vcf.gz
```

---

## Outputs summary

### SV
- `*.precise.vcf.gz` (per-sample)
- `sample_merged_sort.vcf.gz` (svimmer merged)
- `graphtyper.all.vcf.gz`, `graphtyper.PASS.vcf.gz`, `graphtyper.PASS.filtered.vcf.gz`

### SNP/INDEL
- `*.sorted.bam`, `*.sorted.markdup.bam`
- `*.g.vcf.gz`
- `cohort.g.vcf.gz`, `cohort.vcf.gz`
- `cohort.snp.vcf.gz`, `cohort.indel.vcf.gz`
- `cohort.snp.hardfilter.vcf.gz`, `cohort.indel.hardfilter.vcf.gz`
- `cohort.snp.PASS.vcf.gz`, `cohort.indel.PASS.vcf.gz`
### We keep graphtyper.PASS.filtered.vcf.gz and cohort.snp.PASS.vcf.gz files. These two files serve as the foundational dataset for all subsequent filtering and population-level analyses. Refers to variants.vcf.

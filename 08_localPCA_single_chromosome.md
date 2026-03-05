# 08 localPCA (lostruct) for a Single Chromosome (Beagle-imputed SNP VCF)

This document describes a **single-chromosome** localPCA workflow using the `lostruct` R package.
The input is a **Beagle-imputed SNP VCF** (bgzip-compressed), e.g. `beagle.Chr01.edit.vcf.gz`.

Main steps:

1. (Optional) Convert/Filter VCF using **PLINK** (MAF, missingness)
2. Run **localPCA** with `lostruct` using sliding windows along the chromosome
3. Compute window-to-window distances and perform **MDS** (`cmdscale`)
4. Save outputs (`.rds`) and plot MDS1 along genomic position

---

## 0. Input data

- **Single-chromosome** Beagle-imputed SNP VCF (bgzip compressed)  
  Example: `beagle.Chr01.edit.vcf.gz`
- The suffix `.edit.vcf.gz` indicates the VCF has been adjusted to your preferred formatting
  (e.g., chromosome naming, header cleanup, etc.).

Recommended (optional): ensure the VCF is indexed

```bash
tabix -p vcf beagle.Chr01.edit.vcf.gz
```

---

## 1. SNP filtering with PLINK (single chromosome)

### 1.1 Single-file example

```bash
plink --vcf beagle.Chr01.edit.vcf.gz \
  --maf 0.05 \
  --geno 0.5 \
  --recode \
  --transpose \
  --out snp_filter \
  --allow-extra-chr \
  --keep-allele-order
```

Parameter notes:

- `--vcf`: input VCF (Beagle-imputed SNPs)
- `--maf 0.05`: remove SNPs with MAF < 0.05
- `--geno 0.5`: remove SNPs with missing rate > 0.5  
  (for imputed data missingness is usually low; this is relatively permissive)
- `--recode --transpose`: output transposed PLINK files (`.tped/.tfam`)
- `--allow-extra-chr`: allow non-standard chromosome names
- `--keep-allele-order`: keep allele order (avoid PLINK flipping A1/A2 order)

---

## 2. (Optional) Batch PLINK filtering for multiple `*.edit.vcf.gz`

If you have multiple Beagle-imputed VCFs (e.g. multiple chromosomes) in one directory:

```bash
for i in *.edit.vcf.gz
do
  m=${i/.edit.vcf.gz/}
  plink --vcf "$i" \
    --maf 0.05 \
    --geno 0.5 \
    --recode \
    --transpose \
    --out "$m" \
    --allow-extra-chr \
    --keep-allele-order
done
```

> The rest of this document focuses on **one chromosome** (e.g., chr1).

---

## 3. localPCA with `lostruct` (single chromosome)

### 3.1 Load R packages

```r
library(lostruct)
library(data.table)
```

### 3.2 Window the VCF along the chromosome (by physical position)

```r
# size=1e4 means 10,000 bp windows (10 kb)
# type='bp' means windows are defined by physical distance (bp), not by SNP count

snps <- vcf_windower(
  "beagle.Chr01.edit.vcf.gz",
  size = 1e4,
  type = "bp"
)
```

### 3.3 PCA per window (keep first 3 PCs)

```r
pcs <- eigen_windows(snps, k = 3, mc.cores = 64)
```

- `k=3`: compute the first 3 PCs per window
- `mc.cores=64`: number of CPU cores for parallel processing (adjust to your server)

### 3.4 Compute distances between windows (based on PCs)

```r
pcdist <- pc_dist(pcs, npc = 3, mc.cores = 64)
```

- `npc=3`: use the first 3 PCs to compute window-to-window distances

---

## 4. MDS (`cmdscale`) and explained variance

The distance matrix may contain `NA` values (e.g., problematic windows with too few SNPs).
We remove those windows before running MDS.

```r
nas <- is.na(pcdist[, 1])

fit <- cmdscale(pcdist[!nas, !nas], eig = TRUE, k = 3)
```

### 4.1 Explained variance for MDS axes

```r
 eigenvalues <- fit$eig
 explained_variance_total <- sum(eigenvalues)

 explained_variance_MDS1 <- eigenvalues[1] / explained_variance_total
 explained_variance_MDS2 <- eigenvalues[2] / explained_variance_total
 explained_variance_MDS3 <- eigenvalues[3] / explained_variance_total

 print(paste("MDS1 explained variance:", explained_variance_MDS1))
 print(paste("MDS2 explained variance:", explained_variance_MDS2))
 print(paste("MDS3 explained variance:", explained_variance_MDS3))
```

Notes:

- `cmdscale(..., eig=TRUE)` returns eigenvalues.
- If the distance matrix is not perfectly Euclidean, some eigenvalues can be negative.
  If explained variance looks odd, consider summing only positive eigenvalues or reviewing the distance construction.

---

## 5. Save outputs (`.rds`)

`fit$points` are the MDS coordinates for each window.

```r
all <- fit$points

out1 <- paste0("chr1", ".all.rds")
out2 <- paste0("chr1", ".pcdist.rds")
out3 <- paste0("chr1", ".pcs.rds")

saveRDS(all, file = out1)
saveRDS(pcdist, file = out2)
saveRDS(pcs, file = out3)
```

Output meaning:

- `chr1.all.rds`: MDS coordinates (one row per window)
- `chr1.pcdist.rds`: window-to-window distance matrix
- `chr1.pcs.rds`: PCA results per window

---

## 6. Plot MDS1 along genomic position

This section loads the saved MDS coordinates and plots **MDS1** along the chromosome.
Because the sign of MDS axes is arbitrary, we optionally multiply MDS1 by `-1` for visualization.

```r
library(vegan)

x_all <- readRDS("/local_PCA/chr1/chr1.all.rds")

# Optional: flip the sign of MDS1 (MDS direction is arbitrary)
x_all[, 1] <- x_all[, 1] * -1

x_all <- data.frame(x_all)

# Window index is typically stored in rownames
x_all$index <- rownames(x_all)
x_all$index <- gsub('"', "", x_all$index)
x_all$index <- as.numeric(x_all$index)

# Convert window index to Mb.
# With size=10kb windows, index * 10kb = index * 0.01 Mb.
# Dividing by 100 is equivalent to multiplying by 0.01.
x_all$index <- x_all$index / 100

pdf(file = "chr1.pdf")
plot(
  x_all$index, x_all$X1,
  xlab = "Chromosome position (Mb)",
  ylab = "MDS1",
  main = "localPCA: MDS1 along chromosome"
)
dev.off()
```

---

## 7. (Recommended) Make the workflow reusable (optional)

If you plan to run multiple chromosomes, consider parameterizing:

- `chr` name (e.g., `chr1`, `chr19`)
- input VCF path
- window size (`size`)
- number of PCs (`k`, `npc`)
- number of cores (`mc.cores`)
- output prefix

This prevents manual editing errors across chromosomes.

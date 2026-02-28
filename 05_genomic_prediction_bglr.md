# 5. Genomic Prediction (GBLUP/RKHS) with BGLR in R

This script implements genomic prediction using a kernel model in **BGLR** and evaluates prediction accuracy by cross-validation.

---

## Software
- R (>= 3.6 recommended)
- Packages: `data.table`, `BGLR`, `plyr`

---

## Inputs
- Phenotypes: `traits.txt` (tab-delimited, header required)
- Genotypes (imputed SNP/SV matrix): `variants.raw`

## Outputs
- Cross-validation metric table printed to stdout
- Plot of correlation vs fold
- CSV of average metrics: `GBLUPsvgwas_average3.csv`

---

## Step 1 — Convert VCF to raw

```bash
plink --vcf all.missing_maf.recode.vcf \
  --allow-extra-chr --recodeA \
  --out variants
```

---
## Reproducible R script

```r
# -------------------------------
# 0. Load inputs
# -------------------------------
pheno <- read.table("traits.txt", header = TRUE)

library(data.table)
library(BGLR)
library(plyr)

# Genotype matrix: adjust columns to match your .raw format
imputed_snp <- fread("variants.raw")
imputed_snp <- imputed_snp[, -c(1, 3, 4, 5, 6)]  # remove metadata columns (example)

# -------------------------------
# 1. Build standardized genotype matrix and kernel
# -------------------------------
x <- imputed_snp[, -1]                 # remove ID column if present
X <- scale(x, center = TRUE, scale = TRUE)

Y <- read.table("../traits.txt", header = TRUE)
y <- Y[, 1]                            # first trait
n <- length(y)

G <- tcrossprod(X) / ncol(X)           # genomic relationship / kernel

# -------------------------------
# 2. Cross-validation settings
# -------------------------------
nIter  <- 5000
burnIn <- 1000

folds <- 5
set.seed(2024)

# Note: original code mixes folds=5 with "1:10"; here keep 5-fold CV
sets <- rep(1:folds, length.out = n)
sets <- sets[order(runif(n))]

# Metrics containers
COR.CV  <- rep(NA, folds)
MAE.CV  <- rep(NA, folds)
MSE.CV  <- rep(NA, folds)
RMSE.CV <- rep(NA, folds)

yHatCV <- rep(NA, n)

# -------------------------------
# 3. Model fitting per fold
# -------------------------------
for (fold in 1:folds) {
  yNa <- y
  whichNa <- which(sets == fold)
  yNa[whichNa] <- NA

  ETA <- list(list(K = G, model = "RKHS"))

  fm <- BGLR(
    y = yNa,
    ETA = ETA,
    nIter = nIter,
    burnIn = burnIn,
    saveAt = "RKHS_"
  )

  yHatCV[whichNa] <- fm$yHat[fm$whichNa]

  COR.CV[fold]  <- cor(fm$yHat[fm$whichNa], y[whichNa], use = "complete.obs")
  MAE.CV[fold]  <- mean(abs(fm$yHat[fm$whichNa] - y[whichNa]))
  MSE.CV[fold]  <- mean((fm$yHat[fm$whichNa] - y[whichNa])^2)
  RMSE.CV[fold] <- sqrt(MSE.CV[fold])
}

average_results <- data.frame(
  Fold = 1:folds,
  COR  = COR.CV,
  MAE  = MAE.CV,
  MSE  = MSE.CV,
  RMSE = RMSE.CV
)

print(average_results)

plot(average_results$Fold, average_results$COR, type = "b",
     xlab = "Fold", ylab = "Correlation", main = "5-fold CV Accuracy")

write.csv(average_results, "GBLUPsvgwas_average3.csv", row.names = FALSE)
```

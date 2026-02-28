# 02. Population Genetic Analyses

Command-line recipes for population genetic analyses using VCF-based datasets.

## Software
- VCFtools
- PLINK
- ADMIXTURE
- GCTA
- XP-CLR (python implementation)
## Notes
For selection signal analysis, phylogenetic tree construction, and PCA analysis, the file all.missing_maf.recode.vcf (which filters out deletion rates and minor allele frequencies) is used. For structure analysis, the filtered ld file is used because the computation does not depend on ld, significantly improving processing speed. In this paper, we used SNP and SV for admixture and PCA phylogenetic analysis, and FST, π, and XP-CLR were analyzed using quality-controlled SNP data. For details, please see Materials and Methods.

```bash
vcftools --vcf variants.vcf \
  --maf  0.05 \
  --max-missing  0.8 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode --recode-INFO-all \
  --out all.missing_maf
```

## Inputs
- Multi-sample VCF: `all.missing_maf.recode.vcf`
- Sample lists (one sample ID per line): `q3.txt`, `q1.txt`, `q2.txt`, etc.

---

## 2.1 Nucleotide diversity (π), windowed

```bash
vcftools --vcf all.missing_maf.recode.vcf \
  --keep q1.txt \
  --out q1 \
  --window-pi 10000 \
  --window-pi-step 5000
```

---

## 2.2 FST (Weir & Cockerham), windowed

```bash
vcftools --vcf all.missing_maf.recode.vcf \
  --fst-window-size 10000 \
  --fst-window-step 5000 \
  --weir-fst-pop q1.txt \
  --weir-fst-pop q2.txt \
  --keep q1.txt \
  --keep q2.txt \
  --out q1q2
```

---

## 2.3 XP-CLR scan (example Chr01)

```bash
python xpclr \
  --out Chr01.q1q2 \
  --format vcf \
  --input Chr01.vcf \
  --samplesA q1.txt \
  --samplesB q2.txt \
  --chr Chr01 \
  --ld 0.95 \
  --phased \
  --maxsnps 1000 \
  --size 10000 \
  --step 5000
```

---

## 2.4 ADMIXTURE

Convert VCF to BED:

```bash
plink --vcf all.missing_maf.recode.vcf \
  --make-bed --out all \
  --allow-extra-chr \
  --keep-allele-order \
  --set-missing-var-ids @:#
```

Run ADMIXTURE for K=2-4:

```bash
seq 2 4 | awk '{print "admixture --cv -j2 all.bed "$1" 1>admix."$1".log 2>&1"}' > admixture.sh
bash admixture.sh
```

---

## 2.5 PCA (GCTA)

```bash
plink --vcf all.missing_maf.recode.vcf \
  --recode vcf-iid \
  --allow-extra-chr --const-fid \
  --threads 5 \
  --make-bed --out all

gcta64 --bfile all --make-grm --make-grm-alg 0 -autosome --out kinship_yang
gcta64 --grm kinship_yang --pca 20 --out PCA_gcta
```

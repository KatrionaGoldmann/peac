---
title: "Exclude HW snps"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Load in the Data

```{r cars}
library(snpStats)
library(VariantAnnotation)
```


```{r, echo=FALSE}
fl = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/")
fl = fl[grepl("_4PCA.vcf.gz", fl)]
fl = fl[! grepl(".tbi", fl)]
```

## Find the snps which pass the Hardy Weinberg Value


```{}
vcf = readVcf(fl)
res <- genotypeToSnpMatrix(vcf)
gt = res$genotypes

```

```{r}
out = col.summary(gt)
```



```{r}
hardy <- 10^-6      # HWE cut-off
HWEuse <- with(out, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(gt)-sum(HWEuse),"SNPs will be removed due to high HWE.\n") 
cat(sum(HWEuse),"SNPs will be kept.\n") 
out2 = out[HWEuse, ]
````




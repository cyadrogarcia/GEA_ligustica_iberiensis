# Genotype-Environment Association and Outlier Detection Pipeline

This repository contains three scripts used to identify candidate loci under selection or associated with environmental gradients in *Apis mellifera*. Each script performs complementary genome scan analyses using different methods:

1. **SamBada**
2. **RDA (Redundancy Analysis)**
3. **pcadapt**

Input files for each method can be prepared using PLINK

For PCADAPT, regular *bed/*bim/*fam files can be used

For SamBada, *ped, *map in A/C/G/T alleles format can be produced using the following flags --recode --alleleACGT

For RDA, *raw files should be produced using --recodeA. PCA to be used can be produced with --pca

For SamBada and RDA, environmental information for samples should be in the same order as genotype files

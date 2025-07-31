# ============================================================================
# Script: rda_adaptation_analysis.R
# Description: This script performs Redundancy Analysis (RDA) to detect SNPs
#              potentially associated with environmental gradients (local adaptation).
# Author: Carlos A. Yadró Garcia
# Date: 2025-07-31
# Partially adapted from:
# https://github.com/Capblancq/RDA-landscape-genomics
# ============================================================================
# Requirements:
# - R packages: psych, vegan, adegenet, data.table, readxl, robust, qvalue,
#               ggplot2, ggpubr, dplyr
# - Input files:
#     - 'biocl2.csv': Environmental variables
#     - 'iberiensis_maf_geno.eigenvec': PCA of genotypes
#     - 'iberiensis_maf_geno.raw': Genotype matrix from PLINK
#     - 'rdadapt.R': RDAadapt helper function
# ============================================================================

# --- Load libraries ---
library(psych)
library(vegan)
library(adegenet)
library(data.table)
library(readxl)
library(robust)
library(qvalue)
library(ggplot2)
library(ggpubr)
library(dplyr)

# --- Step 1: Load environmental data and PCA coordinates ---
env <- read.csv2("biocl2.csv")
pca <- read.table("iberiensis_maf_geno.eigenvec")
env$PCA1 <- pca$V3
env$PCA2 <- pca$V4
env$PCA3 <- pca$V5
env_red <- env[, 2:ncol(env)]  # remove first column if it's redundant

# --- Step 2: Load and process genotype matrix from PLINK ---
raw_file <- "iberiensis_maf_geno.raw"
matriz <- as.matrix(fread(raw_file, header = TRUE, sep = " "))

# Impute missing values using the most common genotype per SNP
sum(is.na(matriz))
mat.imp <- apply(matriz, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(mat.imp))

# Assign row names
row.names(mat.imp) <- env$ID

# Keep only genotype columns (remove metadata)
mat.red <- mat.imp[, -c(1:6)]

# Convert to numeric matrix
mat.red.num <- matrix(as.numeric(mat.red), ncol = ncol(mat.red)) #ncol(mat.red) should match the number of snps
row.names(mat.red.num) <- row.names(mat.red)
colnames(mat.red.num) <- colnames(mat.red)

# Check consistency
identical(rownames(mat.red.num), env$ID)

# --- Step 3: Run partial RDA with climate variables (controlling for population structure) ---
pRDAclim <- rda(
  mat.red.num ~ 
    lat + 
    lon +
    bio10 + 
    bio13 + 
    bio17 +
    bio2 +
    bio7 +
    bio8 + 
    srad_08.1 +
    Condition(PCA1 + PCA2 + PCA3),
  data = env_red
)

# --- Step 4: Identify candidate SNPs with RDAadapt ---
source("rdadapt.R")
rdadapt_env <- rdadapt(pRDAclim, 3)
range(rdadapt_env$q.values)

# Define FDR threshold
thres_env <- 0.05
outliers005 <- data.frame(
  Loci = colnames(mat.red.num)[which(rdadapt_env$q.values < thres_env)],
  q.value = rdadapt_env$q.values[which(rdadapt_env$q.values < thres_env)],
  contig = unlist(lapply(strsplit(colnames(mat.red.num)[which(rdadapt_env$q.values < thres_env)], "_"), `[`, 1))
)

# --- Step 5: for the candidate SNPs retrieve environmental correlations from the RDA results ---
snps_outliers <- as.character(outliers005$Loci)
ncand <- length(snps_outliers)

foo <- matrix(nrow = ncand, ncol = 10)
colnames(foo) <- colnames(env_red)[1:10]

for (i in 1:ncand) {
  nam <- snps_outliers[i]
  snp.gen <- mat.red.num[, nam]
  foo[i, ] <- apply(env_red[, 1:10], 2, function(x) cor(x, snp.gen))
}
outliers005_cor <- cbind(outliers005, foo)


# Convert correlation columns to numeric
outliers005_cor[, 4:ncol(outliers005_cor)] <- lapply(outliers005_cor[, 4:ncol(outliers005_cor)], as.numeric)

# Identify the environmental variable with highest correlation for each SNP
# the following isntruction will go througth columns for each line to identify the 
# environmental variable with the highest correlation
# in outliers005_cor the environmental variables are in columns 4-13 depending on 
# the number of environmental variable
# which is why  the instruction  goes from column 4 to 13 and creates columns 14 and 15 to store results
outliers005_cor[, 4:ncol(outliers005_cor)] <- lapply(outliers005_cor[, 4:ncol(outliers005_cor)], function(x) as.numeric(as.character(x)))
str(outliers005_cor)
for (i in 1:length(outliers005_cor$Loci)) {
  bar <- outliers005_cor[i,]
  outliers005_cor[i,14] <- names(which.max(abs(bar[4:13]))) # gives the variable
  outliers005_cor[i,15] <- max(abs(bar[4:13]))              # gives the correlation
}


colnames(outliers005_cor)[14] <- "predictor"
colnames(outliers005_cor)[15] <- "correlation"


# Format locus names (remove base pairs if present)
outliers005_cor <- outliers005_cor %>%
  distinct(Loci, .keep_all = TRUE) %>%
  mutate(Loci = sub("_[ATCG]$", "", Loci)) %>%
  distinct(Loci, .keep_all = TRUE)

# --- Step 6: Export results ---
write.csv2(outliers005_cor, file = "outliers005_cor.csv", row.names = FALSE)
writeLines(outliers005_cor$Loci, "snps_outliers005.txt")
write.csv2(outliers005_cor, file = "snpeff/outliers005_cor.csv", row.names = FALSE)
writeLines(outliers005_cor$Loci, "snpeff/snps_outliers005.txt")



#--- Step 7: Create biplots ---

# get outliers pr contig
outliers <- outliers005_cor[order(outliers005_cor$contig, outliers005_cor$q.value), ]
outliers_rdadapt_env <- outliers$Loci[!duplicated(outliers$contig)]

# prepare data for `ggplot`
locus_scores <- scores(pRDAclim, choices = 1:3, display = "species", scaling = "none")
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names %in% outliers005$Loci] <- "qvalue 0.05"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "qvalue 0.05"))
TAB_loci <- TAB_loci[order(TAB_loci$type), ]

TAB_var <- as.data.frame(scores(pRDAclim, choices = 1:3, display = "bp"))

# create biplot
biplot12 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", size = 0.6) +
  geom_point(data = TAB_loci, aes(x = RDA1 * 20, y = RDA2 * 20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "dodgerblue")) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = rownames(TAB_var)), 
            size = 3.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.background = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(size = rel(.8)),
    strip.text = element_text(size = 11),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# Mostrar biplot
print(biplot12)

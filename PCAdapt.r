# Author: Carlos A. Yadró Garcia
# Description: Pipeline to run PCAdapt analysis, visualize outputs, adjust p-values, identify candidate SNPs, and export results.

# ===============================
# Load required libraries
# ===============================
library(pcadapt)
library(qvalue)
library(devtools)
library(OutFLANK)
library(ggplot2)
library(dplyr)

# ===============================
# Load genotype data in BED format
# ===============================
filename <- read.pcadapt("iberiensis_maf_geno.bed", type = "bed")

# ===============================
# Run PCA with K=25 for screeplot visualization
# ===============================
x <- pcadapt(input = filename, K = 25)

# Screeplot of eigenvalues
plot(x, option = "screeplot", K = 25)
pdf("screeplot.pdf")
plot(x, option = "screeplot", K = 25)
dev.off()

# ===============================
# Plot sample scores
# ===============================
plot(x, option = "scores")
pdf("scores.pdf")
plot(x, option = "scores")
dev.off()

# ===============================
# Plot SNP loadings for first 3 PCs
# ===============================
png("loadings.png")
par(mfrow = c(2, 2))
for (i in 1:3) {
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}
dev.off()

# ===============================
# Run PCA again with K = 3 for association testing
# ===============================
x <- pcadapt(filename, K = 3) # define the K according to the screeplot, usually similar to the number of genetic clusters

# Manhattan plot of p-values
png("manhattan.png")
plot(x)
dev.off()

# ===============================
# Adjust p-values
# ===============================
x$qvalues <- qvalue(x$pvalues)$qvalues
x$padjBH <- p.adjust(x$pvalues, method = "BH")
x$padjBF <- p.adjust(x$pvalues, method = "bonferroni")

# Plot distributions of p-value adjustments
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
hist(x$qvalues, xlab = "qvalues", main = NULL, breaks = 50, col = "orange")
hist(x$padjBH, xlab = "padjBH", main = NULL, breaks = 50, col = "orange")
hist(x$padjBF, xlab = "padjBF", main = NULL, breaks = 50, col = "orange")

pdf("qval_hist.pdf")
hist(x$qvalues, xlab = "qvalues", main = NULL, breaks = 50, col = "orange")
dev.off()

pdf("pBH_hist.pdf")
hist(x$padjBH, xlab = "padjBH", main = NULL, breaks = 50, col = "orange")
dev.off()

pdf("qBF_hist.pdf")
hist(x$padjBF, xlab = "padjBF", main = NULL, breaks = 50, col = "orange")
dev.off()

# ===============================
# Create result data frame with SNP IDs and statistics
# ===============================
snps.list <- read.table("iberiensis_maf_geno.bim")

result_df <- data.frame(
  chr = snps.list$V1,
  pos = snps.list$V4,
  snp_id = snps.list$V2,
  pval = x$pvalues,
  qvalues = x$qvalues,
  padjBH = x$padjBH,
  padjBF = x$padjBF
)

# ===============================
# Extract candidate SNPs with qvalue < 0.05
# ===============================
a <- subset(result_df, qvalues < 0.05)
str(a)
length(a$snp_id)

# Save results to CSV and text files
dir.create("snpeff", showWarnings = FALSE)
write.csv2(a, "results_unprunned_q005.csv")
writeLines(a$snp_id, "snp_candidates_unprunned_q005.txt")




# ===============================
# CREAR ManhattanPLOT
# ===============================
library(qqman)
#create with manhattan() from qqman

result_df <- result_df[!is.na(result_df$qvalues), ]
result_df <- result_df[!is.na(result_df$qvalues) & is.finite(result_df$qvalues) & result_df$qvalues > 0, ]


png("manhattanolor2.png", width = 9, height = 6, units = "in", res = 600)
manhattanolor <-manhattan(result_df, chr="chr", bp="pos", snp="snp_id", p="qvalues",
          col = c(adjustcolor("black", alpha.f = 0.5), 
                  adjustcolor("gray40", alpha.f = 0.5)), 
                  suggestiveline = FALSE, 
                  genomewideline = FALSE, 
          chrlabs = c(1:16), ylab="-log10(q)", alpha = 0.5,
          ylim = c(0, 70))
# Agregar las líneas horizontales
abline(h = -log10(0.05), col = "red", lty = 2, lwd = 1.5)   # Línea roja discontinua en -log10(0.05)
abline(h = -log10(0.01), col = "black", lty = 2, lwd = 1.5) # Línea negra discontinua en -log10(0.01)


dev.off()

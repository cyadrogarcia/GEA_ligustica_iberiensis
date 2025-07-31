# ============================================================================
# Script: run_sambada_analysis.R
# Description: This script prepares genotype and environmental data, runs SamBada,
#              performs genomic inflation correction, and extracts candidate SNPs.
# Author: Carlos A. Yadró Garcia
# Date: 2025-07-31
# ============================================================================
# Requirements:
# - R packages: R.SamBada, biomaRt, gdsfmt, htmlwidgets, qvalue, dplyr
# - Input files:
#     - 'iberiensis_maf_geno.ped': PLINK PED file with genotype data
#     - 'biocl2.csv': Environmental variables in CSV format (semicolon-separated)
# - Output directory must be writable
# ============================================================================

# --- Load libraries ---
library(R.SamBada)
library(biomaRt)
library(gdsfmt)
library(htmlwidgets)
library(qvalue)
library(dplyr)

# --- Download SamBada executable ---
downloadSambada(tempdir())  # Optional if already installed
getwd()  # Check current working directory

# --- Step 1: Prepare genotype data ---
prepareGeno(
  pedFile = "iberiensis_maf_geno.ped",
  outputFile = "iberiensis_maf_geno.csv",
  saveGDS = TRUE,
  mafThresh = 0.05,
  missingnessThresh = 0.1,
  interactiveChecks = FALSE
)

# --- Step 2: Prepare environmental data (initial preprocessing) ---
prepareEnv(
  envFile = "biocl2.csv",
  outputFile = "biocl2_sambada.csv",
  maxCorr = 1,
  idName = "ID",
  separator = ";",
  genoFile = "iberiensis_maf_geno.gds",
  numPc = 0.5,
  mafThresh = 0.05,
  missingnessThresh = 0.01,
  ldThresh = NULL,
  numPop = 3,
  clustMethod = "kmeans",
  y = "lat",
  x = "lon",
  interactiveChecks = FALSE,
  verbose = TRUE,
  locationProj = 4326
)

# --- Step 3: Merge environmental and population structure data ---
# As SAMBADA remove highly correlated variables during prepareEnv we are gonna combine the 
# population structure data produced by prepareEnv with the environmental variables that we 
# want to keep (we have previously removed correlated variables)
a <- read.csv2("biocl2.csv", sep = ";")
b <- read.csv2("biocl2_sambada.csv", sep = " ")

new_df <- data.frame(
  ID = a$ID,
  lat = a$lat,
  lon = a$lon,
  bio_10 = a$bio10,
  bio_13 = a$bio13,
  bio_17 = a$bio17,
  bio_2 = a$bio2,
  bio_7 = a$bio7,
  bio_8 = a$bio8,
  bio_9 = a$bio9,
  srad_08 = a$srad_08.1,
  pop1 = b$pop1,
  pop2 = b$pop2
)

# Save as space-separated file
write.table(new_df, file = "biocl3_sambada.csv", sep = " ", row.names = FALSE, quote = FALSE)

# --- Step 4: Run SamBada (parallel implementation) ---
sambadaParallel(
  genoFile = "iberiensis_maf_geno.csv",
  envFile = "biocl3_sambada.csv",
  idGeno = 'ID',
  idEnv = 'ID',
  dimMax = 1,
  saveType = 'END ALL',
  outputFile = 'sambada_exp2',
  populationVar = 'LAST'
)

# --- Step 5: Prepare and explore results ---
prep <- prepareOutput(
  sambadaname = 'sambada_exp2',
  dimMax = 1,
  interactiveChecks = FALSE,
  gdsFile = 'iberiensis_maf_geno.gds',
  popStr = TRUE
)

# --- Step 6: Genomic inflation correction ---
# Based on François, O., Martins, H., Caye, K. & Schoville, S. D. 
# Controlling false discoveries in genome scans for selection. 
# Molecular Ecology 25(454–469), 1365–1294X (2016).

# We can use the gif value directly obtained, but in some times 
# adjustments are needed. See the following optional
z_values <- qnorm(1 - prep$sambadaOutput$pvalueG)
gif <- median(z_values^2) / qchisq(.5, df = 1)
adj.pv <- pchisq(z_values^2 / gif, df = 1, lower.tail = FALSE)

# Optional: Explore histograms for a range of GIF values
# once  we have obtained a gif value from the previous instructions
# we can define a range of gif values to be tested and evaluate the different
#histograms. After that we can decide a gif value to be defined
gif_vals <- seq(0, 15, by = 1)
par(mfrow = c(5, 5), mar = c(4, 4, 2, 1))
for (g in gif_vals) {
  hist(pchisq(z_values^2 / g, df = 1, lower.tail = FALSE),
       col = "green", main = paste("GIF =", g),
       xlab = "Adjusted p-values", xlim = c(0, 1))
}
adj.pv <- pchisq(z_values^2 / gif, df = 1, lower.tail = FALSE)

# Final p- and q-value calibration
adj_q <- qvalue(adj.pv)$qvalues
prep$sambadaOutput$p.adj.G <- adj.pv
prep$sambadaOutput$q.adj.G <- adj_q

# Remove population structure variables
prep$sambadaOutput <- prep$sambadaOutput[!prep$sambadaOutput$Env_1 %in% c("pop1", "pop2"), ]

# --- Step 7: Extract SNP candidates under selection ---

# 7.1: FDR threshold = 0.05
outliers005_cal <- prep$sambadaOutput %>%
  filter(q.adj.G < 0.05) %>%
  group_by(snp) %>%
  slice_min(q.adj.G, with_ties = FALSE) %>%
  ungroup()

write.csv2(outliers005_cal, file = "outliers005_cal.csv", row.names = FALSE)
writeLines(as.character(unique(outliers005_cal$snp)), con = "outliers005_cal.txt")

# ============================================================================
# Step 8: Create Manhattan plots for candidate SNPs using ggplot2
# ============================================================================

# Filter the Sambada output (optional: filter only environmental variables of interest)
test_filtered <- prep$sambadaOutput

# Compute cumulative chromosome positions for plotting
don2 <- test_filtered %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(test_filtered, ., by = c("chr" = "chr")) %>%
  arrange(chr, pos) %>%
  mutate(poscum = pos + tot)

# Compute chromosome center positions for axis labeling
axisdf <- don2 %>%
  group_by(chr) %>%
  summarize(center = (max(poscum) + min(poscum)) / 2)

# Plot Manhattan with facets for each environmental variable
o <- ggplot(don2, aes(x = poscum, y = -log10(qvalueG))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.6, size = 1) +
  scale_color_manual(values = rep(c("black", "gray"), 22)) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    text = element_text(size = 14)
  ) +
  facet_wrap(~ Env_1) +
  geom_hline(yintercept = -log10(0.01), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  xlab("Chromosomes")

# Save the plot to file
ggsave("manhattan_plot.png", plot = o, width = 20, height = 10, units = "in", dpi = 600)


# ============================================================================
# End of script
# ============================================================================

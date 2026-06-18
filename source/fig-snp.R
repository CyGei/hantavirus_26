# Pairwise SNP-distance matrix for the sequenced ANDV L-segment genomes.
#
# Method (core-SNP / complete-deletion): sequences arrive from Pathoplexus
# already aligned to the L-segment reference (no realignment needed); the
# comparison is restricted to the core columns called unambiguously (A/C/G/T) in
# every sequence, and the distance is the raw count of differences over those
# sites (ape model = "N"). See snp_matrix() in R/functions.R.

snp_dist <- snp_matrix(dna)
n_core <- attr(snp_dist, "n_core")   # core positions used for the comparison
n_total <- attr(snp_dist, "n_total") # total aligned positions

pheatmap(
  snp_dist,
  color = colorRampPalette(c("#f7fbff", "#21b528"))(50),
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 20,
  treeheight_col = 20,
  border_color = "white",
  fontsize = 11,
  legend = FALSE,
  main = ""
)

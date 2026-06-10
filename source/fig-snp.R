# =============================================================================
# Pairwise SNP-distance matrix for the sequenced ANDV L-segment genomes
# =============================================================================
# Method (core-SNP / complete-deletion — the standard for a pairwise SNP
# matrix among a small set of closely related consensus genomes):
#
#   1. Alignment. Sequences are retrieved from Pathoplexus already aligned to
#      the ANDV L-segment reference (the alignedNucleotideSequences endpoint),
#      so every genome sits in one fixed coordinate frame and contains no
#      indels. No de-novo realignment (MAFFT/MUSCLE) is therefore needed.
#
#   2. Masking. Consensus bases called below the sequencing-depth threshold are
#      reported as N. We restrict the comparison to the "core" alignment
#      columns that are unambiguously called (A/C/G/T) in *every* sequence.
#      This simultaneously (a) trims the ragged, low-coverage genome ends and
#      any internal low-depth blocks, and (b) guarantees that every pairwise
#      distance is measured over an identical set of sites, so all cells of the
#      matrix are directly comparable (unlike pairwise deletion, where each
#      pair would be scored over a different denominator).
#
#   3. Distance. The SNP distance is the raw count of nucleotide differences
#      between each pair across those core sites (ape's model = "N"). Because
#      the core alignment has no missing data, the pairwise.deletion setting is
#      immaterial here.
# -----------------------------------------------------------------------------

# Aligned sequences -> uppercase character matrix (rows = samples, cols = sites)
aln <- as.character(as.matrix(dna))
aln[] <- toupper(aln)

# Core mask: alignment columns called (A/C/G/T) in *all* sequences
core    <- apply(aln, 2, function(col) all(col %in% c("A", "C", "G", "T")))
n_total <- ncol(aln)   # total aligned positions
n_core  <- sum(core)   # positions retained for the SNP comparison

# Core alignment -> pairwise SNP counts
dna_core <- ape::as.DNAbin(aln[, core, drop = FALSE])
snp_dist <- ape::dist.dna(dna_core, model = "N", pairwise.deletion = FALSE) |>
  as.matrix()

# Order rows/columns by case (WHO) id for a readable matrix
ord      <- order(suppressWarnings(as.integer(rownames(snp_dist))))
snp_dist <- snp_dist[ord, ord]

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

# =============================================================================
# Per-sequence quality summary for the sequenced ANDV L-segment genomes
# =============================================================================
# Reported *before* the SNP-distance calculation so the reader can judge how
# much of each genome was confidently called. Consensus bases below the
# sequencing-depth threshold are reported as N (ambiguous); "completeness" is
# the fraction of aligned positions called A/C/G/T, and the "called span" gives
# the first-to-last called position, exposing ragged low-coverage ends.
# -----------------------------------------------------------------------------

# Aligned sequences -> uppercase character matrix (rows = samples, cols = sites)
aln <- as.character(as.matrix(dna))
aln[] <- toupper(aln)
bases <- c("A", "C", "G", "T")

seq_summary <- map_dfr(rownames(aln), function(id) {
  s <- aln[id, ]
  called_pos <- which(s %in% bases)
  tibble(
    who_id = id,
    length = length(s),
    called = length(called_pos),
    ambiguous = length(s) - length(called_pos),
    completeness = length(called_pos) / length(s),
    span_lo = if (length(called_pos)) min(called_pos) else NA_integer_,
    span_hi = if (length(called_pos)) max(called_pos) else NA_integer_
  )
})

seq_tbl <- seq_summary |>
  left_join(distinct(linelist, who_id, accession_id), by = "who_id") |>
  arrange(suppressWarnings(as.integer(who_id))) |>
  transmute(
    Case = who_id,
    Accession = str_trim(accession_id),
    `Called bases` = format(called, big.mark = ","),
    Ambiguous = format(ambiguous, big.mark = ","),
    Completeness = percent(completeness, accuracy = 0.01),
    `Called span` = str_c(
      format(span_lo, big.mark = ","),
      "–",
      format(span_hi, big.mark = ",")
    )
  )

kable(seq_tbl, align = "llrrrl", booktabs = TRUE)

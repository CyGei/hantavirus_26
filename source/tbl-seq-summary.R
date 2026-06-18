# Per-sequence quality summary for the sequenced ANDV L-segment genomes.
# Reported before the SNP matrix so the reader can judge how much of each genome
# was confidently called. Completeness = fraction of aligned positions called
# A/C/G/T; "called span" = first-to-last called position (exposes ragged ends).

aln <- alignment_matrix(dna) # rows = samples, cols = sites
bases <- c("A", "C", "G", "T")

seq_summary <- map(rownames(aln), \(id) {
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
}) |>
  list_rbind()

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

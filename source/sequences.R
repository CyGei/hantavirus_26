# Aligned ANDV L-segment sequences from Pathoplexus, matched to cases.
# fetch_alignment() caches the download (R/functions.R), so the API is queried
# only when the alignment is absent; subsequent renders run offline.

acc_ids <- linelist |>
  filter(!is.na(accession_id)) |>
  pull(accession_id) |>
  str_trim()

fasta <- here::here(
  "Hondius_hantavirus_h2026", "data", "sequences", "matched_andv_s_alignment.fasta"
)

fetch_alignment(acc_ids, fasta, segment = "L")
dna <- read_alignment(fasta, linelist)

seq_cases <- names(dna) # WHO ids with an associated L-segment sequence

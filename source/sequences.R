acc_ids <- linelist |>
  filter(!is.na(accession_id)) |>
  pull(accession_id) |>
  str_trim()

fasta <- "Hondius_hantavirus_h2026/data/sequences/matched_andv_s_alignment.fasta"

request(
  "https://lapis.pathoplexus.org/andv/sample/alignedNucleotideSequences/L"
) |>
  req_body_json(list(accessionVersion = as.list(acc_ids))) |>
  req_method("POST") |>
  req_perform() |>
  resp_body_raw() |>
  writeBin(fasta)

dna <- read.FASTA(fasta)
names(dna) <- linelist$who_id[match(
  names(dna),
  str_trim(linelist$accession_id)
)]
dna <- dna[!is.na(names(dna))]

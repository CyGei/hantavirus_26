library(httr2)
library(Biostrings)
library(dplyr)

linelist <- read_csv(
  here::here(
    "Hondius_hantavirus_h2026",
    "data",
    "linelist",
    "2026_hantavirus.csv"
  ),
  col_types = cols(
    Gh_ID = col_character(),
    accession_id = col_character(),
    symptom_onset = col_date(),
    outcome_date = col_date(),
    treatment_date = col_date(),
    ship_boarded = col_date(),
    left_ship = col_date(),
    confirmation_date = col_date()
  )
) |>
  filter(status %in% c("confirmed", "probable")) |>
  arrange(symptom_onset) |>
  mutate(
    group = case_when(
      passenger..y.n. == "y" ~ "passenger",
      cruise.crew..y.n. == "y" ~ "crew",
      TRUE ~ NA_character_
    ) |>
      factor(levels = c("passenger", "crew")),
    dates = as.integer(symptom_onset - min(symptom_onset, na.rm = TRUE)),
  )


# 1. Clean and isolate your IDs from your line list dataframe
target_ids <- linelist %>%
  filter(!is.na(accession_id)) %>%
  pull(accession_id) %>%
  unique()

# 2. Define the S-segment alignment endpoint
endpoint_url <- "https://lapis.pathoplexus.org/andv/sample/alignedNucleotideSequences/S"

# 3. Construct and execute the POST request
req <- request(endpoint_url) %>%
  req_body_json(list(accessionVersion = target_ids)) %>%
  req_method("POST")

response <- req_perform(req)

# 4. Save the returned FASTA file locally
writeBin(resp_body_raw(response), "matched_andv_s_alignment.fasta")

writeBin(
  resp_body_raw(response),
  here::here(
    "Hondius_hantavirus_h2026",
    "data",
    "sequences",
    "matched_andv_s_alignment.fasta"
  )
)


# read using ape
library(ape)
dna <- read.FASTA(here::here(
  "Hondius_hantavirus_h2026",
  "data",
  "sequences",
  "matched_andv_s_alignment.fasta"
))
dna

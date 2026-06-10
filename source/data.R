linelist <- read_csv(
  "Hondius_hantavirus_h2026/data/linelist/2026_hantavirus.csv",
  col_types = cols(.default = col_character())
) |>
  clean_names() |>
  filter(case_status %in% c("confirmed", "probable")) |>
  mutate(
    across(matches("date"), mdy),
    who_id = gsub('Case ', '', who_case_number),
    date = as.integer(date_onset - min(date_onset, na.rm = TRUE)),
    group = factor(
      if_else(cruise_passenger_guest == "Y", "passenger", "crew"),
      levels = c("passenger", "crew")
    ),
    ctd_start = ship_board_date,
    ctd_end = pmin(
      ship_disembark_date,
      date_isolation,
      date_death,
      na.rm = TRUE
    )
  ) |>
  select(
    who_id,
    accession_id,
    group,
    date,
    date_onset,
    date_death,
    ctd_start,
    ctd_end,
    contact_setting,
    age,
    gender
  ) |>
  drop_na(date_onset) |>
  arrange(date_onset)

ref_date <- min(linelist$date_onset, na.rm = TRUE)

ctd_timed <- linelist |>
  transmute(
    id = who_id,
    place = contact_setting,
    start = as.integer(ctd_start - ref_date),
    end = as.integer(ctd_end - ref_date),
  )

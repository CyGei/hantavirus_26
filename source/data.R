# Line list and time-resolved shipboard contact windows.
# Logic lives in read_linelist() / build_ctd() (R/functions.R).

linelist <- read_linelist(
  here::here("Hondius_hantavirus_h2026", "data", "linelist", "2026_hantavirus.csv")
)

ref_date <- min(linelist$date_onset, na.rm = TRUE) # day 0 of the integer timeline
ctd_timed <- build_ctd(linelist, ref_date)

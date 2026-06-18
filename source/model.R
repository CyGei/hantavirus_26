# Main reconstruction: symptom onset + shipboard contact windows + L-segment DNA.
# Reconstruction logic lives in R/functions.R; this script wires the inputs and
# exposes the posterior objects used by the figures, tables, and manuscript text.

# Generation-time (w) and incubation-period (f) distributions, assumed equal:
# the UKHSA review reports closely similar incubation (~21 d) and serial
# intervals (~22-23 d) for human-to-human ANDV, so a single discretised gamma
# (mean 22 d, SD 7 d; Martinez 2020 Epuyen, n = 34) is used for both.
si <- serial_interval(mu = 22, cv = 7 / 22)

o2_data <- outbreaker_data(
  ids = linelist$who_id,
  dates = linelist$date,
  dna = dna,
  ctd_timed = ctd_timed,
  w_dens = si,
  f_dens = si,
  p_trans = list(matrix(1, 1, 1))
)

config <- outbreaker_config()

o2_chains <- run_chains(o2_data, config, label = "outbreaker2 chains")
o2 <- bind_chains(o2_chains)
o2_id <- o2ools::identify(o2, ids = linelist$who_id)
entropy_vals <- o2ools::get_entropy(o2_id)
consensus_tree <- o2ools::get_consensus(o2_id)

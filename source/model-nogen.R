# Sensitivity analysis: identical model with the sequence data removed.
# Reuses the serial interval (si) and configuration (config) from source/model.R
# and the same run_chains() helper; the only change is dropping the DNA argument.

o2_data_nogen <- outbreaker_data(
  ids = linelist$who_id,
  dates = linelist$date,
  ctd_timed = ctd_timed,
  w_dens = si,
  f_dens = si,
  p_trans = list(matrix(1, 1, 1))
)

o2_chains_nogen <- run_chains(o2_data_nogen, config, label = "outbreaker2 chains (no sequences)")
o2_nogen <- bind_chains(o2_chains_nogen)
o2_id_nogen <- o2ools::identify(o2_nogen, ids = linelist$who_id)
entropy_vals_nogen <- o2ools::get_entropy(o2_id_nogen)
consensus_tree_nogen <- o2ools::get_consensus(o2_id_nogen)

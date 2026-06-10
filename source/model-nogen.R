# Appendix sensitivity analysis: identical model with the genetic data removed.
# Reuses the serial interval (si) and configuration (config) from source/model.R;
# the only change is dropping the DNA component from the outbreaker data.

o2_data_nogen <- outbreaker_data(
  ids = linelist$who_id,
  dates = linelist$date,
  ctd_timed = ctd_timed,
  w_dens = si,
  p_trans = list(matrix(1, 1, 1))
)

plan(multisession, workers = 4)
o2_chains_nogen <- c(123L, 456L, 789L, 1011L) |>
  future_map(
    \(s) {
      set.seed(s)
      capture.output(
        chain <- outbreaker(data = o2_data_nogen, config = config),
        type = "output"
      )
      filter(chain, step > 10000)
    },
    .options = furrr_options(seed = NULL)
  ) |>
  time_pipe("outbreaker2 chains (no genetics)")
plan(sequential)

o2_nogen <- bind_rows(o2_chains_nogen, .id = "chain")
class(o2_nogen) <- c("outbreaker_chains", class(o2_nogen))
o2_id_nogen <- o2ools::identify(o2_nogen, ids = linelist$who_id)
entropy_vals_nogen <- o2ools::get_entropy(o2_id_nogen)
consensus_tree_nogen <- o2ools::get_consensus(o2_id_nogen)

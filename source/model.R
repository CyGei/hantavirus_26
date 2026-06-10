si <-
  gamma_mucv2shapescale(mu = 18, cv = 6 / 18) |>
  (\(params) {
    distcrete("gamma", 1, shape = params$shape, scale = params$scale)$d(1:50)
  })()

o2_data <- outbreaker_data(
  ids = linelist$who_id,
  dates = linelist$date,
  dna = dna,
  ctd_timed = ctd_timed,
  w_dens = si,
  p_trans = list(matrix(1, 1, 1))
)

config <- create_config(
  n_iter = 1e5,
  sample_every = 50,

  # unobserved intermediates
  max_kappa = 3,
  move_kappa = TRUE,

  # case reporting rate
  init_pi = 0.95,
  move_pi = TRUE,
  prior_pi = c(10, 1),

  # contact reporting coverage
  init_eps = 1,
  move_eps = FALSE,

  # place persitence
  init_tau = 1,
  move_tau = FALSE,

  # importations
  find_import = TRUE,
  n_iter_import = 1e4,

  # mutation rate
  move_mu = FALSE,
  init_mu = 1.7e-5 # ~3.5e-4 subs/site/year × 18/365 days/generation
)


plan(multisession, workers = 4)
o2_chains <- c(123L, 456L, 789L, 1011L) |>
  future_map(
    \(s) {
      set.seed(s)
      capture.output(
        chain <- outbreaker(data = o2_data, config = config),
        type = "output"
      )
      filter(chain, step > 10000)
    },
    .options = furrr_options(seed = NULL)
  ) |>
  time_pipe("outbreaker2 chains")
plan(sequential)

o2 <- bind_rows(o2_chains, .id = "chain")
class(o2) <- c("outbreaker_chains", class(o2))
o2_id <- o2ools::identify(o2, ids = linelist$who_id)
entropy_vals <- o2ools::get_entropy(o2_id)
consensus_tree <- o2ools::get_consensus(o2_id)

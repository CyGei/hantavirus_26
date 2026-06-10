# Sampled parameters of the outbreaker2 model:
#   - pi (scalar, move_pi = TRUE)
#   - t_inf_i (per-case infection times, disctete integers; move_t_inf defaults to TRUE)
# Fixed parameters (mu, eps, tau) are omitted.
diag_cols <- c("pi", stringr::str_subset(names(o2), "^t_inf_\\d+$"))

mcmc_chains <- o2 |>
  group_split(chain) |>
  map(\(df) coda::mcmc(as.matrix(df[diag_cols]))) |>
  coda::mcmc.list()

gd <- coda::gelman.diag(mcmc_chains, autoburnin = FALSE, multivariate = FALSE)
ess <- coda::effectiveSize(mcmc_chains)

n_chains <- n_distinct(o2$chain)
n_samples <- nrow(o2)

diag_tbl <- tibble(
  param = names(ess),
  group = if_else(param == "pi", "pi", "t_inf"),
  rhat = gd$psrf[, "Point est."],
  ess = as.numeric(ess)
) |>
  group_by(group) |>
  summarise(
    n_params = n(),
    `R-hat` = max(rhat),
    ESS = round(min(ess)),
    .groups = "drop"
  ) |>
  arrange(match(group, c("pi", "t_inf"))) |>
  mutate(
    Parameter = if_else(
      group == "pi",
      "Case reporting rate (pi)",
      str_c("Infection times (t_inf, N = ", n_params, ")")
    ),
    `R-hat` = sprintf("%.3f", `R-hat`),
    ESS = format(ESS, big.mark = ",")
  ) |>
  select(Parameter, `R-hat`, ESS)

kable(diag_tbl, align = "lrr", booktabs = TRUE)

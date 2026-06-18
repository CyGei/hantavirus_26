# Cross-chain convergence check: PERMANOVA on the per-chain posterior forests.
# A non-significant result indicates the independent chains explored the same
# distribution of transmission trees. Uses to_forest() from R/functions.R.

o2_forests <- map(o2_chains, \(chain) {
  o2ools::identify(chain, ids = linelist$who_id) |> to_forest()
})

# Seed for reproducibility: permanova_test subsamples trees and permutes labels,
# so an unseeded run gives a p-value that drifts across renders.
set.seed(123)
p_val <- do.call(mixtree::permanova_test, map(o2_forests, sample, size = 200))

p_val |>
  as.data.frame() |>
  tibble::rownames_to_column("Source") |>
  filter(Source == "Model") |>
  transmute(
    df = Df,
    `F` = sprintf("%.2f", `F`),
    `R²` = sprintf("%.3f", R2),
    `p-value` = sprintf("%.3f", `Pr(>F)`),
    Permutations = 999
  ) |>
  kable(align = "rrrrr", booktabs = TRUE)

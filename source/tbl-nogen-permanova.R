# PERMANOVA comparison of the posterior tree distributions with vs without
# sequence data. The two groups are the two pooled posterior forests; a
# significant result would indicate the sequence data shifted the distribution of
# inferred transmission trees. Uses to_forest() from R/functions.R.

forest_gen <- to_forest(o2_id)
forest_nogen <- to_forest(o2_id_nogen)

# Seed for reproducibility: the test subsamples trees and permutes labels, so an
# unseeded run gives a p-value that drifts across renders (it sits near 0.05).
set.seed(123)
p_val_nogen <- mixtree::permanova_test(
  with_sequences = sample(forest_gen, 200),
  without_sequences = sample(forest_nogen, 200)
)

p_val_nogen |>
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

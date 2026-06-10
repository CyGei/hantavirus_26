# PERMANOVA comparison of the posterior tree distributions with vs without
# genetic data, reusing the get_trees() -> mixtree::permanova_test() pattern from
# source/tbl-permanova.R. Here the two groups are the two pooled posterior
# forests; a significant result indicates the genetic data shifted the
# distribution of inferred transmission trees.

to_forest <- function(o2_identified) {
  o2ools::get_trees(o2_identified) |>
    purrr::map(~ filter(.x, !is.na(from)))
}

forest_gen <- to_forest(o2_id)
forest_nogen <- to_forest(o2_id_nogen)

# Seed for reproducibility: the test subsamples trees and permutes labels, so an
# unseeded run gives a p-value that drifts across renders (it sits near 0.05).
set.seed(123)
p_val_nogen <- mixtree::permanova_test(
  with_genetics = sample(forest_gen, 200),
  without_genetics = sample(forest_nogen, 200)
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

# convert each chain into identified chains, then transmission trees and remove the introduction

o2_forests <- purrr::map(.x = o2_chains, .f = \(chain) {
  o2ools::identify(chain, ids = linelist$who_id) |>
    o2ools::get_trees() |>
    purrr::map(~ filter(.x, !is.na(from)))
})

p_val <- do.call(
  mixtree::permanova_test,
  lapply(o2_forests, sample, size = 200)
)

p_val |>
  as.data.frame() |>
  tibble::rownames_to_column("Source") |>
  filter(Source == "Model") |>
  transmute(
    `df` = Df,
    `F` = sprintf("%.2f", `F`),
    `R²` = sprintf("%.3f", R2),
    `p-value` = sprintf("%.3f", `Pr(>F)`),
    `Permutations` = 999
  ) |>
  kable(align = "rrrrr", booktabs = TRUE)

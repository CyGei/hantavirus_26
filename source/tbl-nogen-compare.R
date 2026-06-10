# Per-case comparison of the consensus reconstruction with vs without genetic
# data. Reuses get_consensus() (infector + support) and get_entropy() (ancestry
# uncertainty) outputs already produced for both models. Cases marked * carry an
# L-segment sequence; cases marked † have a different consensus infector when
# genetic data are dropped.

seq_cases <- names(dna)
case_order <- linelist |> arrange(date_onset) |> pull(who_id)

tidy_consensus <- function(consensus, entropy) {
  ent <- tibble(
    to = gsub("alpha_", "", names(entropy)),
    entropy = as.numeric(entropy)
  )
  consensus |>
    transmute(
      to = as.character(to),
      infector = if_else(is.na(from), "Import", as.character(from)),
      support = frequency
    ) |>
    left_join(ent, by = "to")
}

cmp <- tidy_consensus(consensus_tree, entropy_vals) |>
  rename(inf_g = infector, sup_g = support, ent_g = entropy) |>
  left_join(
    tidy_consensus(consensus_tree_nogen, entropy_vals_nogen) |>
      rename(inf_n = infector, sup_n = support, ent_n = entropy),
    by = "to"
  ) |>
  mutate(case = factor(to, levels = case_order), changed = inf_g != inf_n) |>
  arrange(case)

# Summary scalars reused by the appendix narrative.
n_changed <- sum(cmp$changed)
changed_cases <- cmp$to[cmp$changed]
mean_ent_gen <- mean(cmp$ent_g)
mean_ent_nogen <- mean(cmp$ent_n)
mean_ent_seq_gen <- mean(cmp$ent_g[cmp$to %in% seq_cases])
mean_ent_seq_nogen <- mean(cmp$ent_n[cmp$to %in% seq_cases])

cmp |>
  transmute(
    Case = str_c(to, if_else(to %in% seq_cases, "*", ""), if_else(changed, "†", "")),
    Infector = inf_g,
    Support = sprintf("%.2f", sup_g),
    Entropy = sprintf("%.2f", ent_g),
    `Infector ` = inf_n,
    `Support ` = sprintf("%.2f", sup_n),
    `Entropy ` = sprintf("%.2f", ent_n)
  ) |>
  kbl(align = "lcccccc", booktabs = TRUE, escape = FALSE) |>
  add_header_above(c(" " = 1, "With genetics" = 3, "Without genetics" = 3)) |>
  kable_styling(latex_options = "hold_position")

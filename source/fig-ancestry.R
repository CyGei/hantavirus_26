# Ancestry support matrix: for each infectee, how often each infector was
# sampled across the posterior. See ancestry_support() in R/functions.R.
case_order <- linelist |> arrange(date_onset) |> pull(who_id)
from_levels <- c("Import", case_order)

alpha_post <- ancestry_support(o2_id, case_order, from_levels)

consensus_marks <- consensus_tree |>
  mutate(
    from = if_else(is.na(from), "Import", as.character(from)),
    to = factor(to, levels = case_order),
    from = factor(from, levels = from_levels)
  )

ggplot(alpha_post, aes(x = to, y = from)) +
  geom_tile(aes(fill = prob), colour = "grey90", linewidth = 0.2) +
  geom_tile(
    data = consensus_marks,
    fill = NA,
    colour = "#c0392b",
    linewidth = 0.7
  ) +
  scale_fill_gradient(
    "Posterior probability\n",
    low = "white",
    high = "black",
    limits = c(0, 1)
  ) +
  labs(x = "infectee", y = "infector") +
  scale_x_discrete(labels = label_sequenced(seq_cases)) +
  scale_y_discrete(labels = label_sequenced(seq_cases)) +
  coord_fixed() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

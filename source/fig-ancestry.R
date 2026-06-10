case_order <- linelist |> arrange(date_onset) |> pull(who_id)
from_levels <- c("Import", case_order)

alpha_post <- o2_id |>
  select(starts_with("alpha_")) |>
  pivot_longer(
    everything(),
    names_to = "to",
    values_to = "from",
    names_prefix = "alpha_"
  ) |>
  mutate(from = if_else(is.na(from), "Import", as.character(from))) |>
  count(to, from) |>
  group_by(to) |>
  mutate(prob = n / sum(n)) |>
  ungroup() |>
  complete(to = case_order, from = from_levels, fill = list(prob = 0)) |>
  mutate(
    to = factor(to, levels = case_order),
    from = factor(from, levels = from_levels)
  )

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
  coord_fixed() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

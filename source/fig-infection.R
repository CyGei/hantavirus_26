o2 |>
  select(starts_with("t_inf_")) |>
  pivot_longer(
    everything(),
    names_to = "idx",
    names_pattern = "t_inf_(\\d+)",
    names_transform = list(idx = as.integer),
    values_to = "t_inf"
  ) |>
  mutate(
    case = fct_reorder(linelist$who_id[idx], t_inf, mean, .desc = TRUE),
    grp = linelist$group[idx],
    date = ref_date + t_inf
  ) |>
  ggplot(aes(x = date, y = case, fill = grp)) +
  geom_density_ridges(alpha = 0.75, colour = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = group_pal,
    na.value = "grey70",
    labels = group_labels,
    name = NULL
  ) +
  scale_x_date(
    date_breaks = "7 days",
    date_labels = "%d %b",
    expand = expansion(add = 3)
  ) +
  labs(x = "Inferred date of infection", y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

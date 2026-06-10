state_pal <- c(
  "Imported" = "#e6550d",
  "Direct\n(kappa = 1)" = "#525252",
  "1 unobserved\n(kappa = 2)" = "#817f7f",
  "2 unobserved\n(kappa = 3)" = "#c4c3c3"
)

case_order <- linelist |> arrange(desc(date_onset)) |> pull(who_id)

o2 |>
  select(matches("^(alpha|kappa)_\\d+$")) |>
  pivot_longer(
    everything(),
    names_to = c(".value", "idx"),
    names_pattern = "(alpha|kappa)_(\\d+)",
    names_transform = list(idx = as.integer)
  ) |>
  mutate(
    case = factor(linelist$who_id[idx], levels = case_order),
    state = factor(
      case_when(
        is.na(alpha) ~ "Imported",
        kappa == 1 ~ "Direct\n(kappa = 1)",
        kappa == 2 ~ "1 unobserved\n(kappa = 2)",
        .default = "2 unobserved\n(kappa = 3)"
      ),
      levels = names(state_pal)
    )
  ) |>
  count(case, state) |>
  group_by(case) |>
  mutate(prob = n / sum(n)) |>
  ungroup() |>
  complete(case, state, fill = list(n = 0, prob = 0)) |>
  ggplot(aes(x = prob, y = case, fill = state)) +
  geom_col() +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = state_pal, name = NULL) +
  labs(x = "Posterior frequency", y = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

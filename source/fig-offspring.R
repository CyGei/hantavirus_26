get_Ri(o2_id) |>
  mutate(step = row_number()) |>
  pivot_longer(-step, names_to = "name", values_to = "Ri") |>
  count(step, Ri) |>
  group_by(step) |>
  mutate(percent = n / sum(n) * 100) |>
  group_by(Ri) |>
  summarise(
    mean = mean(percent),
    lower = quantile(percent, 0.025),
    upper = quantile(percent, 0.975),
    .groups = "drop"
  ) |>
  ggplot(aes(x = Ri, y = mean)) +
  geom_linerange(
    aes(ymin = 0, ymax = mean),
    colour = "grey60",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_bw() +
  labs(
    x = expression('Number of secondary infections ' ~ R[i]),
    y = "% of cases (posterior mean)"
  )

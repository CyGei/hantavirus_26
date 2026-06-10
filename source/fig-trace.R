param_labels <- c(
  post = "Log posterior",
  pi = "Case reporting rate (pi)"
  # mu = "Mutation rate (mu)",
  # eps = "Contact reporting coverage (eps)"
)

o2 |>
  mutate(chain = paste("Chain", chain)) |>
  select(chain, step, all_of(names(param_labels))) |>
  pivot_longer(-c(chain, step), names_to = "param", values_to = "value") |>
  mutate(param = factor(param_labels[param], levels = param_labels)) |>
  ggplot(aes(x = step, y = value, colour = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.8) +
  facet_wrap(~param, scales = "free_y", ncol = 1) +
  scale_colour_brewer(palette = "Set1", name = NULL) +
  labs(x = "MCMC iteration", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f0f0f0"),
    legend.position = "bottom"
  )

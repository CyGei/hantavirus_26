onset_order <- linelist |> arrange(desc(date_onset)) |> pull(who_id)

tibble(
  case = gsub("alpha_", "", names(entropy_vals)),
  entropy = entropy_vals
) |>
  mutate(case = factor(case, levels = onset_order)) |>
  ggplot(aes(x = entropy, y = case)) +
  geom_col() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 1)) +
  scale_y_discrete(labels = label_sequenced(seq_cases)) +
  labs(x = "Shannon entropy", y = "WHO case ID") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

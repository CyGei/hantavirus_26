ctd_plot <- linelist |>
  mutate(
    sex = str_to_upper(str_sub(gender, 1, 1)),
    label = fct_reorder(
      str_glue("{who_id}  ({sex}, {age})"),
      date_onset
    )
  )

events <- ctd_plot |>
  select(label, date_onset, date_death) |>
  pivot_longer(
    c(date_onset, date_death),
    names_to = "event",
    values_to = "date"
  ) |>
  drop_na(date)

ggplot(ctd_plot) +
  geom_segment(
    aes(
      y = label,
      yend = label,
      x = ctd_start,
      xend = ctd_end,
      colour = group
    ),
    linewidth = 4,
    alpha = 0.75,
    lineend = "round"
  ) +
  geom_point(
    data = events,
    aes(y = label, x = date, colour = event, shape = event),
    size = 3,
    fill = "white"
  ) +
  scale_colour_manual(
    values = c(group_pal, date_onset = "grey20", date_death = "#c0392b"),
    labels = c(
      group_labels,
      date_onset = "Symptom onset",
      date_death = "Death"
    ),
    breaks = c("passenger", "crew", "date_onset", "date_death"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(passenger = NA, crew = NA, date_onset = 21, date_death = 13),
    guide = "none"
  ) +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%d %b",
    expand = expansion(add = 3)
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = c(3, 3, NA, NA),
        shape = c(NA, NA, 21, 13),
        size = c(0, 0, 3, 3),
        fill = c(NA, NA, "white", "white")
      )
    )
  )

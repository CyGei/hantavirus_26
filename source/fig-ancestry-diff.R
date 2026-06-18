# Difference version of the ancestry support matrix (cf. fig-ancestry.R):
# for each infector-infectee pair, the change in posterior support when the
# sequence data are included, P(infector | with sequences) - P(infector | without).
# Green cells are pairs the sequence data favour (support higher with sequences);
# purple cells are pairs it disfavours. Outlines mark each model's consensus
# (modal) infector, coloured to match the fill (dark green = with sequences, dark
# purple = without sequences): a full-cell box for the model with sequences and an inset
# box for the model without sequences, so agreement reads as a box-in-box and a
# disagreement as two separate boxes.

case_order <- linelist |> arrange(date_onset) |> pull(who_id)
from_levels <- c("Import", case_order)

support_gen <- ancestry_support(o2_id, case_order, from_levels) |>
  select(to, from, prob)
support_nogen <- ancestry_support(o2_id_nogen, case_order, from_levels) |>
  select(to, from, prob)

alpha_diff <- support_gen |>
  left_join(
    support_nogen,
    by = c("to", "from"),
    suffix = c("_gen", "_nogen")
  ) |>
  mutate(diff = prob_gen - prob_nogen)

# Consensus (modal) infector under each model, placed on the same factor grid.
model_levels <- c("With sequences", "Without sequences")
mark_consensus <- function(consensus, model) {
  consensus |>
    mutate(
      from = if_else(is.na(from), "Import", as.character(from)),
      to = factor(to, levels = case_order),
      from = factor(from, levels = from_levels),
      model = factor(model, levels = model_levels)
    )
}
consensus_gen <- mark_consensus(consensus_tree, "With sequences")
consensus_nogen <- mark_consensus(consensus_tree_nogen, "Without sequences")

# Symmetric, data-driven limits so the diverging scale is centred at zero.
lim <- max(abs(alpha_diff$diff), na.rm = TRUE)
lim <- if (lim > 0) lim else 1

ggplot(alpha_diff, aes(x = to, y = from)) +
  geom_tile(aes(fill = diff), colour = "grey90", linewidth = 0.2) +
  # with-sequences consensus: full-cell outline (dark green, matching the fill)
  geom_tile(
    data = consensus_gen,
    aes(colour = model),
    fill = NA,
    linewidth = 0.8,
    width = 0.94,
    height = 0.94
  ) +
  # without-sequences consensus: inset outline (dark purple, matching the fill)
  geom_tile(
    data = consensus_nogen,
    aes(colour = model),
    fill = NA,
    linewidth = 0.8,
    width = 0.55,
    height = 0.55
  ) +
  scale_fill_gradient2(
    expression(atop(
      Delta ~ "posterior support",
      "(with" - "without sequences)"
    )),
    low = "#9b45a3",
    mid = "white",
    high = "#4dbd05",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  scale_colour_manual(
    "Consensus infector",
    values = c("With sequences" = "#1b5e20", "Without sequences" = "#5e2766")
  ) +
  labs(x = "infectee", y = "infector") +
  scale_x_discrete(labels = label_sequenced(seq_cases)) +
  scale_y_discrete(labels = label_sequenced(seq_cases)) +
  coord_fixed() +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(7, "lines"),
      barheight = unit(0.6, "lines"),
      order = 1
    ),
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(fill = NA, linewidth = 1),
      order = 2
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "top",
    legend.spacing.x = unit(1.2, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

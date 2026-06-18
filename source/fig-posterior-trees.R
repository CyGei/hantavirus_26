# Four most-supported distinct transmission trees from the posterior.
# A "tree" here is a distinct who-infected-whom ancestry vector (alpha); support
# is the proportion of posterior samples carrying that ancestry. Edges are
# labelled with kappa (the number of transmission generations, 1 = direct).

n_trees <- 4

alpha_cols <- grep("^alpha_", names(o2_id), value = TRUE)
kappa_cols <- grep("^kappa_", names(o2_id), value = TRUE)

# Signature of each posterior sample = its full ancestry vector (NA = import).
sig <- o2_id |>
  select(all_of(alpha_cols)) |>
  mutate(across(everything(), \(x) replace_na(as.character(x), "import"))) |>
  unite("sig", everything(), sep = "|") |>
  pull(sig)

# Rank distinct ancestries by how often they were sampled (posterior support).
n_samples <- length(sig)
tree_rank <- tibble(sig = sig) |>
  count(sig, name = "n", sort = TRUE) |>
  slice_head(n = n_trees) |>
  mutate(tree_id = row_number(), support = n / n_samples)

# Long edge list for the selected trees: one row per (tree, infectee), carrying
# the (constant) infector and the modal kappa across that tree's samples.
samples_long <- o2_id |>
  mutate(.sig = sig) |>
  semi_join(tree_rank, by = c(".sig" = "sig")) |>
  select(.sig, all_of(alpha_cols), all_of(kappa_cols)) |>
  pivot_longer(
    cols = -.sig,
    names_to = c(".value", "to"),
    names_pattern = "(alpha|kappa)_(.*)"
  )

modal <- \(x) {
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

edges <- samples_long |>
  filter(!is.na(alpha)) |> # drop imports (no incoming edge)
  group_by(.sig, to) |>
  summarise(from = first(alpha), kappa = modal(kappa), .groups = "drop") |>
  left_join(tree_rank, by = c(".sig" = "sig"))

# One panel in the consensus-tree style: nodes by onset date (x) and tree
# breadth (y) via the shared tree_layout(); edges coloured by infector role.
plot_one_tree <- function(edges_df, panel_title) {
  layout_data <- tree_layout(edges_df |> select(from, to, kappa), linelist)

  ggraph(layout_data) +
    geom_edge_link(
      aes(color = .N()$group[from]),
      edge_width = 0.6,
      arrow = arrow(length = unit(2.5, "mm")),
      end_cap = circle(3, "mm")
    ) +
    geom_node_point(aes(fill = group), shape = 21, colour = "black", size = 7) +
    geom_node_text(
      aes(label = star_sequenced(name, seq_cases)),
      size = 2.6,
      colour = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(
      NULL,
      values = group_pal,
      breaks = names(group_pal),
      labels = group_labels
    ) +
    scale_edge_colour_manual(values = group_pal, guide = "none") +
    scale_x_continuous(breaks = breaks_width(7), labels = \(x) {
      format(as.Date(x, origin = "1970-01-01"), "%b\n%d")
    }) +
    labs(x = "Symptom Onset Date", y = "", title = panel_title) +
    theme_bw(base_size = 11) +
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
}

panels <- tree_rank$tree_id |>
  map(\(i) {
    plot_one_tree(
      filter(edges, tree_id == i),
      sprintf(
        "Tree %d — %.2f%% of posterior",
        i,
        100 * tree_rank$support[tree_rank$tree_id == i]
      )
    )
  })

wrap_plots(panels, ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

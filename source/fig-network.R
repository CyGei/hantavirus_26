draw_key_arrow <- function(data, params, size) {
  grid::segmentsGrob(
    x0 = 0.1,
    x1 = 0.85,
    y0 = 0.5,
    y1 = 0.5,
    arrow = grid::arrow(length = unit(2, "mm"), type = "closed"),
    gp = grid::gpar(
      col = "grey25",
      fill = "grey25",
      lwd = (data$edge_width %||% 0.5) * ggplot2::.pt
    )
  )
}

epi <- make_epicontacts(
  linelist = linelist,
  contacts = consensus_tree |> drop_na(from),
  id = "who_id",
  directed = TRUE
)

g <- epicontacts:::as.igraph.epicontacts(epi) |> as_tbl_graph()

# Reingold-Tilford: designed for rooted trees, guarantees no overlap
roots <- which(igraph::degree(g, mode = "in") == 0)
tree_layout <- igraph::layout_as_tree(g, root = roots)
rownames(tree_layout) <- igraph::V(g)$name

layout_data <- create_layout(g, layout = "kk")
layout_data$x <- as.numeric(layout_data$date_onset)
layout_data$y <- tree_layout[layout_data$name, 1] # tree breadth → y axis

ggraph(layout_data) +
  geom_edge_link(
    aes(
      edge_width = frequency,
      color = .N()$group[from],
      label = sprintf("%.2f", frequency)
    ),
    arrow = arrow(length = unit(2.5, "mm")),
    end_cap = circle(3, "mm"),
    angle_calc = "along",
    label_dodge = unit(2.5, "mm"),
    label_size = 3,
    label_colour = "black",
    key_glyph = draw_key_arrow
  ) +
  geom_node_point(aes(fill = group), shape = 21, colour = "black", size = 8) +
  geom_node_text(
    aes(label = name),
    size = 2.8,
    colour = "white",
    fontface = "bold"
  ) +
  scale_edge_width(
    "Posterior support",
    range = c(0.1, 1),
    breaks = 1
  ) +
  scale_fill_manual(
    NULL,
    values = group_pal,
    breaks = names(group_pal),
    labels = group_labels
  ) +
  scale_edge_colour_manual(values = group_pal, guide = "none") +
  guides(
    fill = guide_legend(override.aes = list(size = 4), order = 1),
    edge_width = guide_legend(order = 2),
    colour = "none",
    edge_colour = "none",
    edge_alpha = "none"
  ) +
  scale_x_continuous(breaks = breaks_width(3), labels = \(x) {
    format(as.Date(x, origin = "1970-01-01"), "%b\n%d")
  }) +
  labs(x = "Symptom Onset Date", y = "") +
  theme_bw(base_size = 13) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  )

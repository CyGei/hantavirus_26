library(here)
library(tidyverse)
library(outbreaker2)
library(ape)
library(epitrix)
library(distcrete)
library(epicontacts)

# ------------------------------------
#           Data
# ------------------------------------

# --------- Linelist ---------
linelist <- read_csv(
  here::here(
    "Hondius_hantavirus_h2026",
    "data",
    "linelist",
    "2026_hantavirus.csv"
  ),
  col_types = cols(
    Gh_ID = col_character(),
    symptom_onset = col_date(),
    outcome_date = col_date(),
    treatment_date = col_date(),
    ship_boarded = col_date(),
    left_ship = col_date(),
    confirmation_date = col_date()
  )
) |>
  filter(status %in% c("confirmed", "probable")) |>
  arrange(symptom_onset) |>
  mutate(
    group = case_when(
      passenger..y.n. == "y" ~ "passenger",
      cruise.crew..y.n. == "y" ~ "crew",
      TRUE ~ NA_character_
    ) |>
      factor(levels = c("passenger", "crew")),
    dates = as.integer(symptom_onset - min(symptom_onset, na.rm = TRUE)),
  )

group_pal <- c(passenger = "#1f77b4", crew = "#d62728")


# --------- DNA ---------
dna_raw <- ape::read.FASTA(
  here::here(
    "Hondius_hantavirus_h2026",
    "data",
    "sequences",
    "ANDV-Switzerland-Hu-3337-2026.fasta.gz"
  )
)

# --------- CTD ---------
ctd_timed <- linelist |>
  filter(!is.na(ship_boarded), !is.na(left_ship)) |>
  transmute(
    id = Gh_ID,
    place = "MV_Hondius",
    start = as.integer(ship_boarded - min(symptom_onset, na.rm = TRUE)),
    end = as.integer(left_ship - min(symptom_onset, na.rm = TRUE))
  )
ctd_timed
# ------------------------------------
#           outbreaker2
# ------------------------------------
gt_moments <- epitrix::gamma_mucv2shapescale(mu = 18, cv = 6 / 18)
incub_moments <- epitrix::gamma_mucv2shapescale(mu = 18, cv = 6 / 18)

o2_data <- outbreaker_data(
  ids = linelist$Gh_ID,
  dates = linelist$dates,
  #ctd_timed = ctd_timed,
  w_dens = distcrete::distcrete(
    "gamma",
    interval = 1,
    shape = gt_moments$shape,
    scale = gt_moments$scale
  )$d(1:50),
  f_dens = distcrete::distcrete(
    "gamma",
    interval = 1,
    shape = incub_moments$shape,
    scale = incub_moments$scale
  )$d(1:50)
)

config <- create_config(
  n_iter = 1e4,
  sample_every = 50,
  # Allow unobserved intermediate cases between observed ones
  max_kappa = 5,
  move_kappa = TRUE,
  init_pi = 0.9,
  prior_pi = c(1, 1),
  move_pi = TRUE,
  find_import = TRUE
)

set.seed(42)
o2 <- outbreaker(
  data = o2_data,
  config = config
)
# ------------------------------------
#           Results
# ------------------------------------
library(o2ools)
o2 <- o2 |> filter(step > 500)
plot(o2, type = "alpha")
o2_id <- identify(o2, ids = linelist$Gh_ID)
o2_id
entropy <- get_entropy(o2_id)

barplot(
  entropy,
  horiz = TRUE,
  las = 1
)

consensus_tree <- get_consensus(o2_id)
head(consensus_tree)

epi <- make_epicontacts(
  linelist = linelist,
  contacts = subset(consensus_tree, !is.na(from)),
  directed = TRUE
)
plot(epi)


# Requires the `timeline` branch of epicontacts:
#   remotes::install_github("reconhub/epicontacts@timeline")
vis_temporal_interactive(
  epi,
  x_axis = "symptom_onset",
  # network_shape = "rectangle",
  # axis_type = "double",
  date_labels = "%b %d",
  n_breaks = 12,
  thin = FALSE,
  label = "id",
  node_shape = "group",
  shapes = c("passenger" = "user", "crew" = "anchor"), # https://fontawesome.com/v4/icons/
  node_color = "group",
  edge_color = "frequency",
  edge_width = "frequency",
  edge_label = "frequency",
  edge_arrow = "to",
  col_pal = colorRampPalette(group_pal),
  edge_col_pal = colorRampPalette(c("grey80", "grey20"))
)


library(igraph)
library(tidygraph)
library(ggraph)

g <- epicontacts:::as.igraph.epicontacts(epi) |>
  as_tbl_graph()

layout_data <- create_layout(g, layout = 'kk')
layout_data$x <- as.numeric(layout_data$symptom_onset)

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

p <- ggraph(layout_data) +
  geom_edge_link(
    aes(
      edge_width = frequency,
      color = .N()$group[from], # .N() to accesses node data
      label = sprintf("%.2f", frequency)
    ),
    arrow = arrow(length = unit(2.5, 'mm')),
    end_cap = circle(3, 'mm'),
    angle_calc = "along",
    label_dodge = unit(2.5, "mm"),
    label_size = 3,
    label_colour = "black",
    key_glyph = draw_key_arrow
  ) +
  geom_node_point(
    aes(fill = group),
    shape = 21,
    colour = "black",
    size = 5
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3,
    nudge_y = 0.1,
    bg.color = "white",
    bg.r = 0.1
  ) +
  scale_edge_width(
    "Posterior support",
    range = c(0.1, 1),
    breaks = c(0.8, 0.9, 1)
  ) +
  scale_fill_manual(
    NULL,
    values = group_pal,
    breaks = names(group_pal)
  ) +
  scale_edge_colour_manual(
    values = group_pal,
    breaks = names(group_pal),
    guide = "none"
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4), order = 1),
    edge_width = guide_legend(order = 2)
  ) +
  scale_x_continuous(
    breaks = scales::breaks_width(3),
    labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%b %d")
  ) +
  theme_bw() +
  labs(x = "Symptom onset date", y = "") +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  )

print(p)

# Packages, shared constants, and project helpers.
# Grouped by role; only packages actually used downstream are attached.

# Core data manipulation and plotting
library(tidyverse) # dplyr, tidyr, purrr, stringr, ggplot2, forcats, readr, tibble
library(janitor)
library(scales)

# Outbreak reconstruction
library(outbreaker2) # dev version: ctd_timed / p_trans / tau support (see README)
library(o2ools)
library(mixtree)

# Sequences and phylogenetics
library(ape)
library(epitrix)
library(distcrete)

# Transmission-tree graphs
library(epicontacts)
library(tidygraph)
library(ggraph)

# Tables, figures, diagnostics
library(knitr)
library(kableExtra)
library(pheatmap)
library(ggridges)
library(patchwork)
library(coda)

# I/O and parallelism
library(httr2)
library(furrr)
library(future)
library(pipetime)

# Project helpers (pure functions; see R/functions.R)
source(here::here("R", "functions.R"))

# Shared aesthetics
group_pal <- c(passenger = "#1f77b4", crew = "#d62728")
group_labels <- c(passenger = "Passenger", crew = "Crew")

set.seed(123)

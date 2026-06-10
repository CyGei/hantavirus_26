# Running the analysis locally - for debugging and development.
# Call source("run_local.R") in the console.

here::i_am("run_local.R")
setwd(here::here())

FORCE_MODEL <- FALSE

cache_dir <- here::here("_local_cache")
cache_file <- here::here("_local_cache", "model.rds")
dir.create(cache_dir, showWarnings = FALSE)

# 1. Setup, data, sequences (fast) ------------------------------------------
source(here::here("source", "setup.R"))
source(here::here("source", "data.R"))
source(here::here("source", "sequences.R"))

# 2. Model (slow — cached) --------------------------------------------------
# Objects produced by source/model.R that we want available afterwards.
model_objects <- c(
  "si",
  "o2_data",
  "config",
  "o2_chains",
  "o2",
  "o2_id",
  "entropy_vals",
  "consensus_tree"
)

if (!FORCE_MODEL && file.exists(cache_file)) {
  message(
    "Loading cached model from ",
    cache_file,
    " (set FORCE_MODEL <- TRUE to rerun)"
  )
  list2env(readRDS(cache_file), envir = environment())
} else {
  message("Running model...")
  source(here::here("source", "model.R"))
  saveRDS(mget(model_objects), cache_file)
  message("Model cached to ", cache_file)
}

message(
  "\nReady. Globals loaded. Source any figure/table script as needed, e.g.:"
)
message("  source(\"source/fig-network.R\")")

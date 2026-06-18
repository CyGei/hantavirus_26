# =============================================================================
# functions.R — reusable helpers for the MV Hondius ANDV reconstruction
# =============================================================================
# Every function here is pure: it takes explicit inputs and returns an explicit
# value, with no reliance on (or mutation of) global state. The thin scripts in
# source/ wire these together and expose the handful of objects the manuscript's
# inline code references. Packages are attached in source/setup.R; the less
# common ones are namespaced (pkg::fun) to keep each helper's provenance clear.
# =============================================================================


# ---- Line list --------------------------------------------------------------

#' Read and clean the outbreak line list.
#'
#' Keeps confirmed/probable cases that have a symptom-onset date, parses all
#' date columns, and derives the integer day-since-first-onset timeline used by
#' outbreaker2 together with the passenger/crew role and the shipboard exposure
#' window (board date to the earliest of disembarkation, isolation, or death).
#'
#' @param path Path to 2026_hantavirus.csv.
#' @return A tibble, one row per analysable case, arranged by onset date.
read_linelist <- function(path) {
  read_csv(path, col_types = cols(.default = col_character())) |>
    janitor::clean_names() |>
    filter(case_status %in% c("confirmed", "probable")) |>
    mutate(
      across(matches("date"), mdy),
      who_id = gsub("Case ", "", who_case_number),
      date = as.integer(date_onset - min(date_onset, na.rm = TRUE)),
      group = factor(
        if_else(cruise_passenger_guest == "Y", "passenger", "crew"),
        levels = c("passenger", "crew")
      ),
      ctd_start = ship_board_date,
      ctd_end = pmin(ship_disembark_date, date_isolation, date_death, na.rm = TRUE)
    ) |>
    select(
      who_id, accession_id, group, date, date_onset, date_death,
      ctd_start, ctd_end, contact_setting, age, gender
    ) |>
    drop_na(date_onset) |>
    arrange(date_onset)
}

#' Time-resolved contact data (one shared "place" = the ship) for outbreaker2.
#'
#' @param linelist Output of [read_linelist()].
#' @param ref_date Reference date (the earliest symptom onset) defining day 0.
#' @return A tibble with id/place/start/end on the integer-day timeline.
build_ctd <- function(linelist, ref_date) {
  linelist |>
    transmute(
      id = who_id,
      place = contact_setting,
      start = as.integer(ctd_start - ref_date),
      end = as.integer(ctd_end - ref_date)
    )
}


# ---- Sequences --------------------------------------------------------------

#' Fetch the aligned ANDV sequences from Pathoplexus, caching to disk.
#'
#' Downloads only when the FASTA is missing (or `refresh = TRUE`), so that once
#' the alignment has been retrieved, renders are reproducible and work offline.
#' This replaces the previous behaviour of re-POSTing to the API on every render.
#'
#' @param accessions Character vector of accessionVersion ids to request.
#' @param dest Destination FASTA path.
#' @param segment Genome segment endpoint (default "L").
#' @param refresh Force a re-download even if `dest` exists.
#' @return `dest`, invisibly usable as a path.
fetch_alignment <- function(accessions, dest, segment = "L", refresh = FALSE) {
  if (refresh || !file.exists(dest)) {
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    httr2::request(paste0(
      "https://lapis.pathoplexus.org/andv/sample/alignedNucleotideSequences/",
      segment
    )) |>
      httr2::req_body_json(list(accessionVersion = as.list(accessions))) |>
      httr2::req_method("POST") |>
      httr2::req_perform() |>
      httr2::resp_body_raw() |>
      writeBin(dest)
  }
  dest
}

#' Read a cached alignment and relabel sequences by WHO case id.
#'
#' Sequences whose accession does not match a case in the line list (e.g. a
#' malformed or unresolved accession) are dropped.
#'
#' @param path FASTA path.
#' @param linelist Output of [read_linelist()].
#' @return A `DNAbin` list named by `who_id`.
read_alignment <- function(path, linelist) {
  dna <- ape::read.FASTA(path)
  names(dna) <- linelist$who_id[match(names(dna), str_trim(linelist$accession_id))]
  dna[!is.na(names(dna))]
}


# ---- Alignment / SNP distances ----------------------------------------------

#' Aligned `DNAbin` -> uppercase character matrix (rows = samples, cols = sites).
alignment_matrix <- function(dna) {
  aln <- as.character(as.matrix(dna))
  aln[] <- toupper(aln)
  aln
}

#' Core pairwise SNP-distance matrix (complete deletion).
#'
#' Restricts the comparison to the "core" alignment columns called
#' unambiguously (A/C/G/T) in *every* sequence, then counts raw nucleotide
#' differences (ape model = "N"). Restricting to shared, high-confidence sites
#' trims ragged low-coverage ends and guarantees every pairwise distance is
#' measured over an identical set of positions. Rows/columns are ordered by WHO
#' id; `n_core` (sites used) and `n_total` (aligned positions) are attached as
#' attributes for the manuscript's inline reporting.
#'
#' @param dna A `DNAbin` list of aligned sequences.
#' @return A numeric matrix with attributes `n_core` and `n_total`.
snp_matrix <- function(dna) {
  aln <- alignment_matrix(dna)
  core <- apply(aln, 2, function(col) all(col %in% c("A", "C", "G", "T")))
  d <- ape::as.DNAbin(aln[, core, drop = FALSE]) |>
    ape::dist.dna(model = "N", pairwise.deletion = FALSE) |>
    as.matrix()
  ord <- order(suppressWarnings(as.integer(rownames(d))))
  d <- d[ord, ord]
  attr(d, "n_core") <- sum(core)
  attr(d, "n_total") <- ncol(aln)
  d
}


# ---- Transmission model -----------------------------------------------------

#' Discretised gamma serial-interval density over `1:max_days`.
#'
#' @param mu Mean (days).
#' @param cv Coefficient of variation.
#' @param max_days Support upper bound.
serial_interval <- function(mu, cv, max_days = 50) {
  p <- epitrix::gamma_mucv2shapescale(mu = mu, cv = cv)
  distcrete::distcrete("gamma", 1, shape = p$shape, scale = p$scale)$d(seq_len(max_days))
}

#' outbreaker2 configuration shared by the main and no-genetics models.
#'
#' Defined once so the sensitivity analysis is guaranteed to differ from the
#' main analysis only in the data supplied. `init_mu` (2.1e-5 subs/site/
#' generation) is the published ANDV rate of ~3.5e-4 subs/site/year scaled to
#' the 22-day generation time.
outbreaker_config <- function() {
  outbreaker2::create_config(
    n_iter = 1e5,
    sample_every = 50,
    max_kappa = 3,        # up to 2 unobserved intermediates
    move_kappa = TRUE,
    init_pi = 0.95,       # case reporting rate
    move_pi = TRUE,
    prior_pi = c(10, 1),
    init_eps = 1,         # contact reporting coverage (fixed)
    move_eps = FALSE,
    init_tau = 1,         # place persistence (fixed)
    move_tau = FALSE,
    find_import = TRUE,   # importation detection
    n_iter_import = 1e4,
    move_mu = FALSE,      # mutation rate fixed (near-zero diversity)
    init_mu = 2.1e-5
  )
}

#' Run independent outbreaker2 chains in parallel; return post-burn-in chains.
#'
#' Each worker sets its own seed, so the run is reproducible; furrr's parallel
#' RNG stream is therefore disabled (`seed = NULL`). The parallel plan is reset
#' on exit. Returns the list of per-chain tibbles (kept separately so chain-wise
#' diagnostics such as the PERMANOVA can use them).
#'
#' @param data An `outbreaker_data` object.
#' @param config An `outbreaker_config` object.
#' @param seeds Integer seeds, one per chain.
#' @param burnin Iterations discarded from the start of each chain.
#' @param workers Parallel workers.
#' @param label Label for the timing message.
run_chains <- function(data, config,
                       seeds = c(123L, 456L, 789L, 1011L),
                       burnin = 10000, workers = 4,
                       label = "outbreaker2 chains") {
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  seeds |>
    furrr::future_map(
      \(s) {
        set.seed(s)
        capture.output(
          chain <- outbreaker2::outbreaker(data = data, config = config),
          type = "output"
        )
        filter(chain, step > burnin)
      },
      .options = furrr::furrr_options(seed = NULL)
    ) |>
    pipetime::time_pipe(label)
}

#' Pool per-chain tibbles into a single `outbreaker_chains` object.
bind_chains <- function(chains) {
  o2 <- bind_rows(chains, .id = "chain")
  class(o2) <- c("outbreaker_chains", class(o2))
  o2
}

#' Posterior forest from an identified chain.
#'
#' Converts an identified chain to its list of transmission trees and removes
#' the introduction edge (`from = NA`) from each, the form `mixtree` expects.
to_forest <- function(o2_identified) {
  o2ools::get_trees(o2_identified) |>
    map(\(tree) filter(tree, !is.na(from)))
}

#' Posterior ancestry-support matrix, in long form.
#'
#' For each infectee (`to`), the proportion of posterior samples (`prob`) that
#' assign each potential infector (`from`); imports (`alpha = NA`) are labelled
#' "Import". Completed over the full `case_order` x `from_levels` grid so missing
#' pairs are explicit zeros, and returned with both columns as ordered factors.
#'
#' @param o2_identified An identified `outbreaker_chains` object.
#' @param case_order Character vector of infectee ids (x-axis order).
#' @param from_levels Character vector of infector levels, e.g. c("Import", case_order).
#' @return A tibble with columns `to`, `from`, `n`, `prob`.
ancestry_support <- function(o2_identified, case_order, from_levels) {
  o2_identified |>
    select(starts_with("alpha_")) |>
    pivot_longer(
      everything(),
      names_to = "to", values_to = "from", names_prefix = "alpha_"
    ) |>
    mutate(from = if_else(is.na(from), "Import", as.character(from))) |>
    count(to, from) |>
    group_by(to) |>
    mutate(prob = n / sum(n)) |>
    ungroup() |>
    complete(to = case_order, from = from_levels, fill = list(prob = 0)) |>
    mutate(
      to = factor(to, levels = case_order),
      from = factor(from, levels = from_levels)
    )
}


# ---- Labels -----------------------------------------------------------------

#' Append "*" to case ids that carry a sequence.
#'
#' Marks sequenced cases on a plot without altering the underlying ids (so joins
#' and factor ordering still use the bare `who_id`). Use directly on a vector of
#' ids — e.g. ggraph node labels — or via [label_sequenced()] for a discrete axis.
#'
#' @param ids Character/factor vector of case ids.
#' @param seq_ids Ids that carry a sequence (e.g. `seq_cases`).
star_sequenced <- function(ids, seq_ids) {
  paste0(as.character(ids), if_else(as.character(ids) %in% seq_ids, "*", ""))
}

#' Labeller factory for `scale_*_discrete(labels = label_sequenced(seq_cases))`.
label_sequenced <- function(seq_ids) {
  function(ids) star_sequenced(ids, seq_ids)
}


# ---- Transmission-tree layout -----------------------------------------------

#' Rooted-tree layout for a transmission graph.
#'
#' Nodes are placed by symptom-onset date (x) and Reingold-Tilford tree breadth
#' (y, which guarantees no edge overlap for a rooted tree). Shared by the
#' consensus-tree figure and the posterior-trees panel so both use an identical
#' layout. `edges` must carry `from`/`to` plus any edge attributes the caller
#' maps in ggraph (e.g. `frequency` or `kappa`).
#'
#' @param edges Edge tibble (from, to, ...).
#' @param linelist Output of [read_linelist()] (supplies node attributes).
#' @return A ggraph `layout_data` data frame.
tree_layout <- function(edges, linelist) {
  epi <- make_epicontacts(
    linelist = linelist, contacts = edges, id = "who_id", directed = TRUE
  )
  g <- epicontacts:::as.igraph.epicontacts(epi) |> as_tbl_graph()

  # Reingold-Tilford layout from the introduction(s) (in-degree 0).
  roots <- which(igraph::degree(g, mode = "in") == 0)
  rt <- igraph::layout_as_tree(g, root = roots)
  rownames(rt) <- igraph::V(g)$name

  layout_data <- create_layout(g, layout = "kk")
  layout_data$x <- as.numeric(layout_data$date_onset) # onset date -> x axis
  layout_data$y <- rt[layout_data$name, 1]            # tree breadth -> y axis
  layout_data
}

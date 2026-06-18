# hantavirus_26

Bayesian inference of transmission chains for the 2026 Andes virus (ANDV) outbreak linked to the MV *Hondius* cruise ship.

## What it does

`manuscript.qmd` renders a self-contained HTML/PDF manuscript that:

1. Reads the case line list from [kraemer-lab/Hondius_hantavirus_h2026](https://github.com/kraemer-lab/Hondius_hantavirus_h2026) (git submodule)
2. Fetches aligned ANDV L-segment sequences from [Pathoplexus](https://pathoplexus.org) (cached after first download)
3. Probabilistically infers who infected whom using **outbreaker2**, summarising the posterior as a consensus tree, ancestry/entropy support, an offspring distribution, and convergence diagnostics

## Project layout

```
manuscript.qmd        Quarto manuscript; each analysis step is a chunk that
                      sources one script from source/ via `#| file:`
R/functions.R         Reusable, side-effect-free helpers (data, sequences,
                      SNP distances, model runner, tree layout)
source/               Thin orchestration scripts (setup, data, sequences,
                      model, model-nogen) and one script per figure/table
run_local.R           Run the pipeline interactively, caching the model
references.bib         Bibliography (CSL: eid.csl)
Hondius_hantavirus_h2026/   Data submodule (line list + sequences)
```

## Clone

```bash
git clone --recurse-submodules https://github.com/CyGei/hantavirus_26
```

> `--recurse-submodules` is required to populate the `Hondius_hantavirus_h2026/` data folder.

## Render

```r
# In R, from the project directory
quarto::quarto_render("manuscript.qmd")
```

Or from the terminal:

```bash
quarto render manuscript.qmd
```

For interactive development, `source("run_local.R")` runs setup → data → sequences → model and caches the (slow) model fit to `_local_cache/`, so figure and table scripts can be sourced individually without refitting.

## Dependencies

R packages are attached in `source/setup.R`. The reconstruction relies on the
development versions of **outbreaker2** (time-resolved contact data: `ctd_timed`,
`p_trans`, `tau`), **o2ools**, and **mixtree**, which are not on CRAN. Pin exact
sources/commits with [`renv`](https://rstudio.github.io/renv/) for full
reproducibility, e.g.:

```r
renv::init()
# install the dev packages, then:
renv::snapshot()   # writes renv.lock
```

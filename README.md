# hantavirus_26

Bayesian inference of transmission chains for the 2026 Andes virus (ANDV) outbreak linked to the MV *Hondius* cruise ship.

## What it does

`manuscript.qmd` renders a self-contained HTML/PDF manuscript that:

1. Reads the case line list from [kraemer-lab/Hondius_hantavirus_h2026](https://github.com/kraemer-lab/Hondius_hantavirus_h2026) (git submodule)
2. Fetches aligned ANDV L-segment sequences from [Pathoplexus](https://pathoplexus.org) (cached after first download)
3. Probabilistically infers who infected whom using **outbreaker2**.

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
`p_trans`, `tau`), **o2ools**, and **mixtree**.

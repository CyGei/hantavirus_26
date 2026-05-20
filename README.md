# hantavirus_26

Bayesian transmission tree reconstruction for the 2026 Andes virus (ANDV) outbreak linked to the MV Hondius cruise ship.

## What it does

`report.qmd` renders a self-contained HTML/PDF report that:

1. Pulls the latest case data from [kraemer-lab/Hondius_hantavirus_h2026](https://github.com/kraemer-lab/Hondius_hantavirus_h2026)
2. Fetches aligned L-segment sequences from [Pathoplexus](https://pathoplexus.org)
3. Runs 4 parallel MCMC chains with **outbreaker2** to infer who infected whom
4. Produces transmission trees, infection date distributions, convergence diagnostics, and import/kappa probabilities

## Clone

```bash
git clone --recurse-submodules https://github.com/CyGei/hantavirus_26
```

> `--recurse-submodules` is required to populate the `Hondius_hantavirus_h2026/` data folder.

## Render

```r
# In R, from the project directory
quarto::quarto_render("report.qmd")
```

Or from the terminal:

```bash
quarto render report.qmd
```

## Dependencies

All R packages are loaded in the `setup` chunk of `report.qmd`. Key packages: `outbreaker2`, `o2ools`, `furrr`, `ggraph`, `gt`, `ggridges`, `httr2`, `coda`.

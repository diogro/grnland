# AGENTS.md

Guidance for coding agents working in this repository. For design and
internals, read `DEVELOPMENT.md`.

## What this is

`grnland` is an R package (R + Rcpp/RcppArmadillo) that simulates populations
whose phenotypes are steady states of gene regulatory networks, evolving on
adaptive landscapes. It generalises the code behind the paper *Developmental
Bias Shapes Evolution on Complex Adaptive Landscapes* (Netz, Melo & Milocco).

## Layout

| Path | Contents |
|---|---|
| `R/model.R` | model constructors: `grn()`, `paper_motif()`, `grn_custom()`, locus parsing |
| `R/develop.R` | `develop()`, `phenotype()`, `sensitivity()`, `m_matrix()`, `genotype_at()` |
| `R/linear.R` | `gp_linear()`: linear genotype-phenotype map |
| `R/landscape.R` | landscape constructors, `landscape_grid()`, `find_peaks()` |
| `R/simulate.R` | `simulate_population()`, selection schemes, stopping rules |
| `R/lande.R` | `lande_separatrices()` |
| `src/hill.cpp` | Hill-network solver and analytical derivatives (the hot path) |
| `tests/testthat/` | tests; `helper-legacy.R` holds closed forms from the original scripts |
| `vignettes/grnland.Rmd` | getting-started vignette (built by R CMD check) |
| `vignettes/articles/` | pkgdown-only article with the paper figures (not built by check) |
| `analysis/` | scripts reproducing the paper figures; not part of the package |
| `docs/` | built pkgdown site, served by GitHub Pages; committed |

## Commands

Run from the repository root. Cap BLAS threads for anything that runs
simulations.

```sh
R CMD INSTALL .                                         # reinstall after any change
Rscript -e 'devtools::test()'                           # tests (~1 min)
Rscript -e 'devtools::document()'                       # after editing roxygen comments
Rscript -e 'Rcpp::compileAttributes()'                  # after changing exported C++ signatures
Rscript -e 'pkgdown::build_site(preview = FALSE)'       # rebuild docs/ (installs the package)
R CMD build . && R CMD check --no-manual grnland_*.tar.gz
```

- Analysis scripts and vignettes use the *installed* package. Reinstall
  before running them, or they will silently use the old code.
- Run analysis scripts from `analysis/` with
  `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`; test them with
  `QUICK=1` first. Full runs take many core-hours; ask before launching them.
- `R CMD check` gives one WARNING on the maintainer's Mac about a missing
  Homebrew gcc search path. It comes from the local toolchain; anything else
  is new.

## Conventions

- Edit roxygen comments in `R/`, never `man/*.Rd` or `NAMESPACE` directly.
- Never edit `R/RcppExports.R` or `src/RcppExports.cpp`; regenerate them.
- `regulation[i, j]` is the effect of gene `j` on gene `i`. Edges are written
  `"a -> b"` (a activates b) and `"a -| b"` (a represses b).
- Locus codes shared by R and C++: 1 = alpha, 2 = gamma, 3 = theta, 4 = hill.
  C++ indices are 0-based.
- Development status codes: 0 ok, 1 no convergence, 2 unstable, 3 non-finite.
- New behaviour gets a test. Derivative code gets a finite-difference check
  (see `test-general.R`); motif results are checked against the closed forms
  in `helper-legacy.R`.
- Keep defaults equal to the settings of the paper's scripts (haploid clonal,
  N = 1000, mu = 0.01, sigma = 0.05, `selection_relative_min()`), unless the
  change is deliberate and documented in the README.
- Simulations take a `seed` and must be reproducible from it.

## Docs

`docs/` is committed and served at https://diogro.github.io/grnland/. After
changing exported functions, the README, the vignette or the article, run
`pkgdown::build_site()` and commit `docs/` together with the change. The
article's figures in `vignettes/articles/figures/` are copied from
`analysis/figures/`; most are from reduced `QUICK=1` runs, and the article
says so.

## Don'ts

- Don't commit `analysis/output/` or `analysis/figures/` (cached results and
  generated figures; both ignored).
- Don't change the random-number call order in `simulate_population()` without
  a reason: it changes every seeded result, including the analysis caches.
- Don't push, or change repository settings, unless asked.

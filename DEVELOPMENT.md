# Developing grnland

This document describes how the package works and how to extend it. For a
user-level introduction, see the README and `vignette("grnland")`.

## Setup

You need R, a C++ compiler, and the packages in `DESCRIPTION`
(`Rcpp`, `RcppArmadillo`, plus `ambient`, `ggplot2`, `knitr`, `rmarkdown` and
`testthat` for the suggested features). For development, also install
`devtools`, `roxygen2` and `pkgdown`.

```sh
git clone git@github.com:diogro/grnland.git
cd grnland
R CMD INSTALL .
Rscript -e 'devtools::test()'
```

## Architecture

A simulation chains four pieces, each replaceable on its own:

```
model  ──►  develop(model, genotypes, x_start)  ──►  phenotypes  ──►  landscape(X)  ──►  fitness
  ▲                                                                                       │
  └──────────── simulate_population(): selection, reproduction, mutation ◄────────────────┘
```

### Models (`R/model.R`, `R/linear.R`)

A model is a list with a class, and these fields: `genes`, `n_genes`,
`traits` (indices of the selected genes), `loci` (a data frame with one row
per evolvable parameter: `name`, `code`/`param`, `i`, `j`) and `x0` (the
default starting expression state).

| Class | Constructor | Development |
|---|---|---|
| `hill_grn` | `grn()`, `paper_motif()` | C++ (`src/hill.cpp`) |
| `custom_grn` | `grn_custom()` | R port of the same solver (`custom_develop_one()`) |
| `linear_gp` | `gp_linear()` | closed form, `X = G B'` |

`grn()` parses edges into a regulation matrix (`parse_regulation()`), expands
`evolvable` into loci (`hill_loci()`), and resolves `traits`.

### Development (`R/develop.R`, `src/hill.cpp`)

`develop()` is an S3 generic returning `list(X, status)`: the full expression
state of every individual and a status code (0 ok, 1 no convergence,
2 unstable, 3 non-finite).

For Hill networks, `hill_develop_cpp()` takes the base parameters, the locus
table, a genotype matrix and starting states, and for each individual:

1. writes the genotype into the network parameters (`apply_genotype()`);
2. integrates the ODE with adaptive Dormand–Prince 5(4) until the right-hand
   side is small (1e-5);
3. polishes the steady state with Newton's method to `tol`. It accepts the
   result only if it is non-negative and stays near the integrated state, so
   it cannot jump to another branch. Otherwise it integrates further, to `tol`;
4. checks stability from the eigenvalues of the Jacobian.

Starting states matter for multistable networks (e.g. motif 5): during
simulations every offspring starts development from its parent's state, so
populations stay on their branch unless mutation pushes them off it.

### Derivatives

`hill_derivs_cpp()` returns the right-hand side `f`, the Jacobian `J` and
`df/dp` for every locus, all analytical. From these,

- `sensitivity()`: `S = dx*/dp = -J^{-1} df/dp` (implicit function theorem);
- `m_matrix()`: `M = S_T Σ S_T'`, where `S_T` holds the rows of the trait
  genes and `Σ = ploidy · μ · diag(σ²)`;
- `genotype_at()`: inverts the map. When the loci are exactly the alphas, the
  steady-state equations are linear in alpha and are solved directly;
  otherwise it runs Newton's method on the loci using `S`.

Custom models use the user's Jacobian if given, otherwise finite differences.

### Landscapes (`R/landscape.R`)

A landscape is a function from an `n × d` phenotype matrix to `n` fitness
values, with class `grn_landscape` and attributes `n_dim` and `info`. The
attribute is `n_dim` and not `dim`, because `structure()` refuses to set
`dim` on a function. `landscape_grid()` precomputes a 2D landscape and
interpolates bilinearly. With Perlin noise this matters: `ambient` has a
large fixed cost per call, so caching makes simulations about three times
faster.

`find_peaks()` scans a grid in 2D, merging maxima closer than `merge` cells.
That stops flat tops from producing duplicate peaks. In other dimensions it
runs multistart `optim()`.

Perlin noise from `ambient` is exactly flat along x within some lattice cells,
on the lines `y · frequency = k + ½`. Expect degenerate critical points there.
`lande_separatrices()` ignores them and points at the domain edge.

### Simulation (`R/simulate.R`)

State: alleles `A` (`N × L × ploidy`; for diploids each copy carries half the
starting value), expression `X` (`N × n_genes`), fitness `W`, and a flag for
individuals whose development failed. Each generation:

1. **Selection.** `selection(W)` turns fitness into sampling weights. The
   default `selection_relative_min()` is the original scripts' `w - min(w) + 0.001`.
2. **Reproduction.** Clonal reproduction resamples individuals. Sexual
   reproduction samples mothers and fathers (from separate halves when
   `sexes = TRUE`) and builds gametes with free recombination.
3. **Mutation.** Draws the number of mutations from `Binomial(N·L·ploidy, μ)`,
   adds `Normal(0, σ)` effects, and clamps the alleles to `[lower, upper]`.
4. **Development.** Only individuals whose genotype changed are developed
   again, starting from their parent's expression state. Failed developments
   are lethal (fitness NA, weight 0) or, with `invalid = "reject"`, the
   mutation is undone.
5. **Recording and stopping.** Every `record_every` generations it stores
   means, the phenotypic covariance and optionally M at the mean genotype.
   Every `stop$every` generations the stopping rule is checked.

The cost is dominated by development of changed individuals and by landscape
evaluation. At the paper's settings a generation takes about 1 ms with a
cached landscape.

### Lande separatrices (`R/lande.R`)

For deterministic evolution `dz/dt = G ∇W`, with `G ∝ M(z)` (or a fixed `G`),
it works as follows:
1. Evaluate the landscape on a grid and take finite-difference gradients and
   Hessians.
2. Find critical points by Newton's method from a lattice of seeds, and keep
   the non-degenerate saddles.
3. Trace the stable manifold of each saddle backwards with RK4.

Because `G` is positive definite, the critical points and their types are the
same as those of the landscape.

## Extending

**A new locus type for Hill networks** (e.g. a basal-expression term):
1. Add the parameter to the dynamics in `rhs()` and `jac()` in `hill.cpp`.
2. Give it a new code in `apply_genotype()` and its derivative in `dfdp()`.
3. Handle the new class name in `hill_loci()`.
4. Add a finite-difference test of `sensitivity()`.

**A new model class:**
1. Write a constructor returning the fields listed under *Models*.
2. Write a `develop.<class>()` method.
3. If `sensitivity()`, `genotype_at()` and `default_genotype()` can't use the
   generic code, add branches for the class. Everything downstream, including
   `m_matrix()` and `simulate_population()`, then works unchanged.

**A new landscape:** wrap any function with `landscape_function(f, dim)`, or
write a constructor that calls `new_landscape()`.

**A new selection scheme or stopping rule:** a selection scheme is a function
from fitness to weights. A stopping rule is a list with `every` and
`fun(X) -> peak index or NULL`. See `stop_at_peaks()`.

## Testing

- `test-motifs.R`: the general solver against the closed forms of the paper's
  motifs (steady states to about 1e-8, and the sensitivities of `Fig2.r`).
- `test-general.R`: edge parsing, loci, finite-difference derivative checks,
  custom vs Hill models, linear maps, `lande_separatrices()`.
- `test-simulate.R`: hill-climbing, reproducibility from a seed, diploid
  sexual runs, `find_peaks()`, `landscape_grid()` accuracy.

Keep tests fast (the whole suite runs in about a minute) and seeded.

## Documentation

- Roxygen comments in `R/` generate `man/` and `NAMESPACE`: run
  `devtools::document()`.
- `vignettes/grnland.Rmd` is a real vignette, built by `R CMD check`; it must
  run in seconds.
- `vignettes/articles/paper-figures.Rmd` is pkgdown-only (excluded in
  `.Rbuildignore`). It shows code without running it and embeds images from
  `vignettes/articles/figures/`, which are copied from `analysis/figures/`.
  When the full figure runs are done, copy the new images over.
- The site: `pkgdown::build_site()` writes `docs/`, which is committed and
  served by GitHub Pages from `main`/`docs`. Equations use MathJax
  (`math-rendering: mathjax` in `_pkgdown.yml`).

## Paper analyses

`analysis/` reproduces the paper figures; `analysis/README.md` lists scripts,
settings and run times. Scripts share `common.R`, cache results in
`analysis/output/` and accept `QUICK=1`, `N_CORES`, `REDO=1`. Always cap BLAS
threads (`OPENBLAS_NUM_THREADS=1` etc.) when running them in parallel.

## Releasing

1. Bump `Version` in `DESCRIPTION`.
2. Run `devtools::document()`, `devtools::test()`, then
   `R CMD build . && R CMD check --no-manual grnland_*.tar.gz`.
3. Rebuild the site with `pkgdown::build_site()` and commit `docs/`.
4. Tag the release (`git tag v0.x.y`) and push the tag.

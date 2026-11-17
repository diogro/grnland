# grnland

Individual-based simulations of populations whose phenotypes are the steady
states of gene regulatory networks (GRNs), evolving on arbitrary adaptive
landscapes. This package generalises the simulation code behind *Developmental
Bias Shapes Evolution on Complex Adaptive Landscapes* (Netz, Melo & Milocco) to
networks with any number of genes and any regulatory structure.

## Installation

```r
# from a local clone
devtools::install("path/to/grnland")
```

Requires a C++ compiler (Rcpp/RcppArmadillo). `ambient` is needed for Perlin
landscapes.

## Models

### Hill-type networks: `grn()`

For genes $i = 1..n$,

$$\dot x_i = \alpha_i \prod_j H_{ij}(x_j) - \gamma_i x_i,$$

where $H_{ij}(x) = x^h/(\theta^h + x^h)$ if gene $j$ activates gene $i$,
$1 - x^h/(\theta^h + x^h)$ if $j$ represses $i$, and 1 if they don't
interact ($\theta = \theta_{ij}$, $h = h_{ij}$). The phenotype is the steady
state that development reaches.

```r
library(grnland)

# edges: "a -> b" activation, "a -| b" repression; self-regulation allowed
m <- grn(c("a -> b", "b -| c", "c -> a", "a -| c"),
         alpha = c(4, 6, 8), theta = 3,
         evolvable = c("alpha", "theta"),   # which parameters are loci
         traits = c("b", "c"))              # which genes are under selection
m
#> <hill_grn> 3 genes, 7 loci, traits: b, c
#> edges: a -> b, a -| c, b -| c, c -> a
#> loci: alpha[a], alpha[b], alpha[c], theta[a,c], theta[b,a], theta[c,a], theta[c,b]
```

`regulation` can also be a matrix where `regulation[i, j]` is the effect of
gene `j` on gene `i` (+1, -1 or 0). `evolvable` accepts whole parameter
classes (`"alpha"`, `"gamma"`, `"theta"`, `"hill"`) or single entries
(`"alpha[b]"`, `"theta[c,a]"`).

The six motifs of the paper are available as `paper_motif(1)` to `paper_motif(6)`.

### Arbitrary dynamics: `grn_custom()`

For regulatory logic that the Hill model cannot express (additive inputs,
basal expression, other kinetics), give the right-hand side of the ODE
directly. This runs in R, so it is slower.

```r
m <- grn_custom(
  function(x, p) c(p[["a1"]] - x[1],
                   p[["a2"]] * (0.2 + x[1] / (5 + x[1])) - x[2]),
  params = c(a1 = 2, a2 = 2), genes = c("x1", "x2"))
```

## Genotype-phenotype map

```r
g <- default_genotype(m)            # or genotype_at(m, phenotype)
phenotype(m, g)                     # steady-state trait values
develop(m, g)                       # all genes plus status (stable / failed)
sensitivity(m, g)                   # J, df/dp and S = dx*/dp = -J^{-1} df/dp
m_matrix(m, g, mu = 0.01, sigma = 0.05, ploidy = 1)   # M = S Sigma S'
```

Development integrates the ODE (Dormand–Prince 5(4)) from a starting state,
polishes the steady state with Newton's method and checks that it is stable.
In networks with several stable states, the starting state (`x_start`) decides
which one is reached; during simulations each offspring starts from its
parent's state.

## Landscapes

A landscape is any function from a phenotype matrix (one row per individual)
to fitness values. There are constructors for common cases:

```r
landscape_perlin(frequency = 0.5, seed = 1245)            # 2D or 3D Perlin noise
landscape_gaussian(rbind(c(4, 2), c(2, 5)), omega = 1)    # Gaussian peaks, any dimension
landscape_smoothmax(l1, l2, beta = 100)                   # smooth maximum of landscapes
landscape_grid(landscape_perlin(0.5, 1245))               # cached 2D grid (~3x faster sims)
find_peaks(land, lower = 0, upper = 15)                   # local optima
```

## Simulation

```r
land <- landscape_grid(landscape_perlin(frequency = 0.5, seed = 1245))
peaks <- find_peaks(land, 0, 15)
sim <- simulate_population(
  paper_motif(1), land, start = c(4, 6),
  N = 1000, generations = 10000, mu = 0.01, sigma = 0.05,
  ploidy = 1, reproduction = "clonal",          # or ploidy = 2, "sexual"
  selection = selection_relative_min(),         # or selection_proportional(), selection_exponential(s)
  stop = stop_at_peaks(peaks), record_m = TRUE, seed = 1)
sim$outcome       # index of the peak reached (NA if none)
sim$trajectory    # means, phenotypic (co)variances P_*, M_*, mean fitness per recorded generation
```

The defaults reproduce the settings of the paper's scripts: haploid clonal
populations of 1000, per-locus mutation rate 0.01, mutational sd 0.05, alleles
truncated at 0, the `w - min(w) + 0.001` selection weights and the "half the
population within ±0.1 of a peak" convergence rule.
`inst/examples/grid_screen.R` shows a full grid screen of starting points.

## Differences from the original analysis scripts

- Steady states are computed numerically for every network, not from closed
  forms. For the six motifs they agree with the closed forms to about 1e-8,
  and the sensitivities agree with the analytical expressions in `Fig2.r`
  (see `tests/testthat/test-motifs.R`).
- Motif 5 (mutual activation): when $\alpha_1\alpha_2 < \theta^2$ the network
  goes to the zero state, which is what the dynamics do; the closed form
  gives negative expression there. The stale-phenotype bug of `Fig5.r` does
  not exist here: phenotypes are always recomputed from the mutated genotype.
- Individuals whose development fails (no stable steady state) get zero
  fitness by default (`invalid = "lethal"`); `invalid = "reject"` discards
  the offending mutation instead.
- Random numbers are drawn differently, so runs reproduce the paper's results
  statistically, not draw for draw.

## Tests

```r
devtools::test()
```

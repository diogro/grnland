# Adaptive landscapes

A landscape is a function taking a matrix of phenotypes (one row per
individual) and returning a vector of fitness values. Any such function
can be used directly; these constructors cover the common cases.

## Usage

``` r
landscape_perlin(
  frequency,
  seed,
  offset = 0,
  transform = function(z) (z + 1)/2,
  dim = 2,
  ...
)

landscape_gaussian(optima, omega = 1, heights = 1, combine = c("sum", "max"))

landscape_smoothmax(..., beta = 100)

landscape_function(f, dim)

landscape_grid(landscape, lower = 0, upper = 15, step = 0.005)
```

## Arguments

- frequency:

  Perlin frequency; higher values give more rugged landscapes.

- seed:

  Perlin seed.

- offset:

  Added to phenotypes before evaluating the noise (scalar or one value
  per dimension), e.g. to shift the landscape off the lattice.

- transform:

  Applied to the raw noise in \\\[-1, 1\]\\; the default rescales to
  \\\[0, 1\]\\.

- dim:

  Number of phenotypic dimensions.

- ...:

  Further arguments to \`ambient::gen_perlin()\`.

- optima:

  Matrix of peak positions, one row per peak.

- omega:

  Width of the peaks: scalar variance, vector of variances, or a
  covariance matrix (shared by all peaks), or a list of matrices.

- heights:

  Peak heights.

- combine:

  \`"sum"\` or \`"max"\` over peaks.

- beta:

  Sharpness; large values approach the maximum.

- f:

  Function taking a phenotype matrix and returning fitness values.

- landscape:

  A 2D landscape to cache.

- lower, upper:

  Grid bounds (scalar or per dimension).

- step:

  Grid spacing.

## Functions

- `landscape_perlin()`: Perlin noise (via the ambient package), in 2 or
  3 dimensions.

- `landscape_gaussian()`: Sum (or maximum) of Gaussian peaks, \\W(x) =
  \sum_k h_k \exp(-\frac12 (x-\theta_k)^\top \Omega^{-1}
  (x-\theta_k))\\.

- `landscape_smoothmax()`: Smooth maximum of several landscapes,
  \\\sum_k w_k e^{\beta w_k} / \sum_k e^{\beta w_k}\\.

- `landscape_function()`: Wrap any function of a phenotype matrix.

- `landscape_grid()`: Precompute a 2D landscape on a grid and
  interpolate bilinearly. Much faster for landscapes that are slow to
  evaluate (e.g. Perlin noise, where each call has a large fixed cost),
  at the price of a small interpolation error. Outside the grid, values
  are clamped to the nearest edge.

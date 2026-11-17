# Develop genotypes into phenotypes

Integrates the regulatory dynamics from a starting expression state
until a steady state is reached. The steady state is checked for
stability.

## Usage

``` r
develop(
  model,
  genotypes,
  x_start = NULL,
  tol = 1e-10,
  tmax = 10000,
  maxsteps = 1e+05,
  check_stability = TRUE
)
```

## Arguments

- model:

  A \`grn\` model.

- genotypes:

  Matrix of genotypic values, one row per individual and one column per
  locus (see \[loci_names()\]), or a single vector.

- x_start:

  Starting expression states (matrix with one row per individual, or one
  vector used for all). Defaults to \`model\$x0\`. When the network has
  several stable states, this determines which is reached.

- tol:

  Relative tolerance on \\\|dx/dt\|\\ at the steady state.

- tmax, maxsteps:

  Limits on developmental time and integration steps.

- check_stability:

  Check that the steady state is stable.

## Value

A list with \`X\`, the full gene-expression matrix, and \`status\` (0 =
stable steady state, 1 = no convergence, 2 = unstable steady state, 3 =
non-finite values).

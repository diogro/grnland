# Separatrices of the Lande-equation flow

Treats evolution of the mean phenotype as the flow \\\dot{\bar z} =
G(\bar z) \nabla W(\bar z)\\, with the G-matrix taken proportional to
the local M-matrix of the genotype-phenotype map, and traces the
separatrices between the basins of attraction of different optima: the
stable manifolds of the saddle points of the landscape. Because G is
positive definite, critical points and their types are those of the
landscape; only the separatrices move. Two dimensions only.

## Usage

``` r
lande_separatrices(
  model,
  landscape,
  lower = 0,
  upper = 12,
  step = 0.01,
  G = NULL,
  seed_step = 0.5,
  ds = 0.01,
  max_steps = 4000
)
```

## Arguments

- model:

  A \`grn\` model with two traits (used for the default \`G\`).

- landscape:

  A 2D landscape function.

- lower, upper:

  Bounds of the region (scalar or per dimension).

- step:

  Grid spacing for the landscape derivatives.

- G:

  Function of a phenotype returning the 2x2 G-matrix. Defaults to the
  model's M-matrix at that phenotype.

- seed_step:

  Spacing of the starting points used to find critical points.

- ds, max_steps:

  Arc-length step and maximum number of steps when tracing.

## Value

Data frame of separatrix paths with columns \`x1\`, \`x2\`, \`saddle\`,
\`branch\` and \`step\`. The saddle points are in \`attr(, "saddles")\`.

## Details

The landscape is evaluated once on a grid; gradient and Hessian come
from finite differences on that grid and are interpolated bilinearly (as
in the original analysis scripts). Only the direction of the flow
matters, so G can be given up to a scalar factor, and \\\nabla \ln W\\
and \\\nabla W\\ give the same separatrices.

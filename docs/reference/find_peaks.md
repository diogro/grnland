# Find local optima of a landscape

In two dimensions, evaluates the landscape on a regular grid and returns
cells that are at least as high as all eight neighbours, dropping the
grid border (as in the original analysis scripts). In other dimensions,
runs \`optim()\` from random starts and merges nearby solutions.

## Usage

``` r
find_peaks(
  landscape,
  lower = 0,
  upper = 15,
  res = 500,
  dim = attr(landscape, "n_dim"),
  n_starts = 200,
  merge = 3
)
```

## Arguments

- landscape:

  A landscape function.

- lower, upper:

  Bounds of the search region (scalar or per dimension).

- res:

  Grid resolution per dimension (2D).

- dim:

  Number of dimensions (taken from the landscape if possible).

- n_starts:

  Number of random starts (other dimensions).

- merge:

  In 2D, maxima closer than this many grid cells are merged, keeping the
  highest (flat tops otherwise give duplicate peaks). 0 keeps all, as in
  the original scripts.

## Value

Data frame with peak coordinates and \`fitness\`.

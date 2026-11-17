# Derivatives of the developmental system at a steady state

Derivatives of the developmental system at a steady state

## Usage

``` r
sensitivity(model, genotype, x = NULL)
```

## Arguments

- model:

  A \`grn\` model.

- genotype:

  Vector of genotypic values.

- x:

  Steady-state expression of all genes. Computed if omitted.

## Value

A list with \`x\`, the Jacobian \`J\` (\\\partial f/\partial x\\),
\`dfdp\` (\\\partial f/\partial\\ loci) and the sensitivity matrix \`S\`
\\= \partial x^\*/\partial\\ loci \\= -J^{-1} \partial f/\partial p\\
(all genes by all loci).

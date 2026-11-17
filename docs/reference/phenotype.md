# Phenotype (trait values) of genotypes

Convenience wrapper around \[develop()\] returning the trait columns
only. Individuals whose development fails get \`NA\`.

## Usage

``` r
phenotype(model, genotypes, x_start = NULL, ...)
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

- ...:

  Passed to \[develop()\].

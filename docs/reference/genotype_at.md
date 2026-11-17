# Genotype that develops into a given phenotype

Analytical for Hill networks where every gene is a trait and the loci
are exactly the production rates \`alpha\`. Otherwise solved by Newton
iteration on the sensitivity matrix, which needs as many loci as traits.

## Usage

``` r
genotype_at(model, phenotype, start = NULL, tol = 1e-10)
```

## Arguments

- model:

  A \`grn\` model.

- phenotype:

  Target trait values.

- start:

  Starting genotype for the Newton iteration.

- tol:

  Tolerance on the phenotype.

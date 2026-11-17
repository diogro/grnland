# Mutational variance-covariance matrix (M-matrix)

\$\$M = S_T \Sigma S_T^\top,\$\$ where \\S_T\\ holds the rows of the
sensitivity matrix belonging to the traits and \\\Sigma\\ is the
per-generation mutational variance of genotypic values, \`ploidy \* mu
\* sigma^2\` per locus. This is the local (linearised) M-matrix at the
given genotype.

## Usage

``` r
m_matrix(
  model,
  genotype = NULL,
  phenotype = NULL,
  mu = 0.01,
  sigma = 0.05,
  ploidy = 1
)
```

## Arguments

- model:

  A \`grn\` model.

- genotype:

  Genotypic values. Either this or \`phenotype\` is needed.

- phenotype:

  Trait values; the genotype is found with \[genotype_at()\].

- mu:

  Per-allele mutation rate per locus (scalar or per locus).

- sigma:

  Standard deviation of mutational effects (scalar or per locus).

- ploidy:

  1 (haploid) or 2 (diploid).

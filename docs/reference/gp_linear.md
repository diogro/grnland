# Linear genotype-phenotype map

Phenotype \\x = B \alpha\\, the linear map used in Fig 1 of the paper
(Wagner 1984). Its M-matrix is constant: \\M = \Sigma B B^\top\\ with
\\\Sigma\\ = \`ploidy \* mu \* sigma^2\`.

## Usage

``` r
gp_linear(B, genes = paste0("x", seq_len(nrow(B))))
```

## Arguments

- B:

  Square matrix mapping genotypic values to traits.

- genes:

  Trait names.

## Value

An object of class \`linear_gp\`, usable wherever a \`grn\` model is.

## Examples

``` r
B2 <- matrix(c(5, 4, 4, 5) / 3, 2)
m <- gp_linear(B2)
m_matrix(m, c(1, 1), mu = 0.01, sigma = 0.05)
#>              x1           x2
#> x1 0.0001138889 0.0001111111
#> x2 0.0001111111 0.0001138889
```

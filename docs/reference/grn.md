# Hill-type gene regulatory network

Builds a network of \`n\` genes with dynamics \$\$dx_i/dt = \alpha_i
\prod_j H\_{ij}(x_j) - \gamma_i x_i,\$\$ where \\H\_{ij}(x) =
x^h/(\theta^h + x^h)\\ if gene j activates gene i, \\1 - x^h/(\theta^h +
x^h)\\ if j represses i, and 1 otherwise (with \\\theta =
\theta\_{ij}\\, \\h = h\_{ij}\\). The phenotype is the steady state
reached by development.

## Usage

``` r
grn(
  regulation,
  genes = NULL,
  alpha = 1,
  gamma = 1,
  theta = 1,
  hill = 1,
  evolvable = "alpha",
  traits = NULL,
  x0 = 1
)
```

## Arguments

- regulation:

  Either an \`n x n\` matrix with entries -1, 0 or 1, where
  \`regulation\[i, j\]\` is the effect of gene \`j\` on gene \`i\` (+1
  activation, -1 repression), or a character vector of edges such as
  \`c("x1 -\> x2", "x2 -\| x1")\`.

- genes:

  Gene names, or the number of genes. Defaults to the row names of
  \`regulation\`, or the genes named in the edges (in order of
  appearance).

- alpha, gamma:

  Production and decay rates (scalar or length \`n\`). \`alpha\` values
  are only defaults when alpha is evolvable.

- theta, hill:

  Regulatory thresholds and Hill coefficients (scalar or \`n x n\`
  matrix, indexed like \`regulation\`).

- evolvable:

  Which parameters are genetically variable loci. Any of \`"alpha"\`,
  \`"gamma"\`, \`"theta"\`, \`"hill"\` (all entries of that kind; for
  theta and hill only existing edges), or single entries such as
  \`"alpha\[2\]"\`, \`"theta\[2,1\]"\`, mixed freely.

- traits:

  Genes whose expression is the phenotype under selection (names or
  indices). Defaults to all genes.

- x0:

  Default initial expression state used when no other starting point is
  given.

## Value

An object of class \`hill_grn\`.

## Examples

``` r
# motif 3 of the paper: x2 represses x1, x1 activates x2
m <- grn(c("x1 -> x2", "x2 -| x1"), theta = 5)
# a three-gene repressilator-like cascade with evolvable thresholds
m3 <- grn(c("a -| b", "b -| c", "c -> a"), evolvable = c("alpha", "theta"))
```

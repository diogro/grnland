# User-defined regulatory network

For regulatory relations that the Hill model in \[grn()\] cannot
express. Development and derivatives run in R, so this is slower than
\[grn()\].

## Usage

``` r
grn_custom(
  rhs,
  params,
  genes,
  evolvable = names(params),
  traits = NULL,
  jacobian = NULL,
  x0 = 1
)
```

## Arguments

- rhs:

  Function \`rhs(x, p)\` returning \\dx/dt\\ for expression vector \`x\`
  and named parameter vector \`p\`.

- params:

  Named numeric vector of all parameters (defaults).

- genes:

  Gene names, or number of genes.

- evolvable:

  Names of parameters in \`params\` that are loci.

- traits:

  Genes under selection (names or indices).

- jacobian:

  Optional function \`jacobian(x, p)\` returning the \\n \times n\\
  matrix \\\partial f/\partial x\\. Numerical if omitted.

- x0:

  Default initial expression state.

## Value

An object of class \`custom_grn\`.

## Examples

``` r
# additive (OR-logic) activation of x2 by x1, which grn() cannot express
m <- grn_custom(
  function(x, p) c(p["a1"] - x[1],
                   p["a2"] * (0.2 + x[1] / (5 + x[1])) - x[2]),
  params = c(a1 = 2, a2 = 2), genes = c("x1", "x2"))
```

# Keys for de-duplicated development

\`dev_key_exact()\` keys on the exact genotype and the exact starting
expression state, so de-duplicated development gives bit-identical
results to developing every individual. \`dev_key_genotype()\` keys on
the genotype only. Each genotype is then developed once from the
starting state of its first occurrence. That is also exact when each
genotype has a single stable state reachable from all the starting
states present (the usual case away from bistable networks), and
otherwise can pick a different steady state. Values are compared exactly
(hexadecimal representation of the doubles).

## Usage

``` r
dev_key_exact(G, X0)

dev_key_genotype(G, X0)
```

## Arguments

- G:

  Genotype matrix, one row per individual.

- X0:

  Starting expression states, one row per individual.

## Value

A character vector with one key per row.

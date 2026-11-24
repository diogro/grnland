# Develop each distinct individual once

Wraps a development function so that individuals with the same key are
developed only once and share the result. This saves time when many new
individuals share a genotype, as with discrete allelic values (few
alleles per locus, recombined in sexual populations). With continuous
mutational effects nearly every genotype is distinct, and it only adds
overhead.

## Usage

``` r
develop_unique(dev, key = dev_key_exact)
```

## Arguments

- dev:

  A function \`dev(G, X0)\` of a genotype matrix and a matrix of
  starting expression states, one row per individual.

- key:

  A function \`key(G, X0)\` returning one string per individual.
  Individuals with equal keys must develop identically. See
  \[dev_key_exact()\] and \[dev_key_genotype()\].

## Value

A function with the same interface as \`dev\`.

## Details

The wrapper is independent of the simulation loop. It applies to any
function \`dev(G, X0)\` that returns a list whose elements have one row
(or one entry) per individual, as \[develop()\] does (\`X\`,
\`status\`). Other kinds of individuals, for example several motifs in
one population, need a \`key\` that tells them apart (e.g. a motif id),
or no de-duplication.

## See also

\[simulate_population()\] (argument \`dedupe\`).

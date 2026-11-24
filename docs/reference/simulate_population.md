# Simulate a population evolving on an adaptive landscape

Individual-based simulation with discrete, non-overlapping generations:
selection (parents sampled with fitness-dependent weights),
reproduction, mutation, then development of new genotypes into
phenotypes.

## Usage

``` r
simulate_population(
  model,
  landscape,
  start = NULL,
  start_genotype = NULL,
  N = 1000,
  generations = 10000,
  mu = 0.01,
  sigma = 0.05,
  ploidy = 1,
  reproduction = c("clonal", "sexual"),
  sexes = TRUE,
  selection = selection_relative_min(),
  lower = 0,
  upper = Inf,
  invalid = c("lethal", "reject"),
  record_every = 10,
  record_m = FALSE,
  record_genotype = FALSE,
  record_third_moments = FALSE,
  dedupe = c("none", "exact", "genotype"),
  stop = NULL,
  seed = NULL,
  dev_args = list()
)
```

## Arguments

- model:

  A \`grn\` model.

- landscape:

  Function from a phenotype matrix to fitness values.

- start:

  Starting phenotype (trait values). The population starts monomorphic
  at the genotype found by \[genotype_at()\].

- start_genotype:

  Alternatively, the starting genotypic values.

- N:

  Population size.

- generations:

  Maximum number of generations.

- mu:

  Per-allele mutation rate per locus.

- sigma:

  Standard deviation of mutational effects (scalar or per locus).

- ploidy:

  1 or 2.

- reproduction:

  \`"clonal"\` or \`"sexual"\`.

- sexes:

  For sexual reproduction, split the population into equal numbers of
  females and males.

- selection:

  A selection scheme, see \[selection\].

- lower, upper:

  Bounds on allelic values.

- invalid:

  What happens to individuals whose development fails (no stable steady
  state): \`"lethal"\` gives them zero weight; \`"reject"\` (clonal
  only) discards the mutation that caused the failure.

- record_every:

  Record population summaries every this many generations (the first and
  last generations are always recorded).

- record_m:

  Also record the M-matrix at the mean genotype.

- record_genotype:

  Also record the mean genotypic value of each locus
  (\`gmean\_\<locus\>\`) and the genotypic covariance matrix of the loci
  (\`G\_\<locus\>\_\<locus\>\`, upper triangle), over individuals that
  develop.

- record_third_moments:

  Also record the third central moments of the genotypic values,
  E\[(g_i - mean_i)(g_j - mean_j)(g_k - mean_k)\] for i \<= j \<= k
  (\`G3\_\<locus\>\_\<locus\>\_\<locus\>\`; divisor n, not n - 1), over
  individuals that develop.

- dedupe:

  Develop each distinct new individual once per generation and share the
  result (\[develop_unique()\]). Only worth it when genotypes repeat
  often, e.g. with discrete allelic values; with continuous mutational
  effects nearly every genotype is unique and it only adds overhead.
  \`"none"\` (default) develops every new individual. \`"exact"\` keys
  on the genotype and the starting expression state and gives the same
  results as \`"none"\`. \`"genotype"\` keys on the genotype only, which
  saves more but can pick a different steady state for genotypes with
  several stable states (see \[dev_key_genotype()\]).

- stop:

  Optional stopping rule, see \[stop_at_peaks()\].

- seed:

  Optional RNG seed.

- dev_args:

  List of extra arguments for \[develop()\], e.g. \`list(newton_radius =
  0)\` to always integrate instead of starting development of new
  individuals with Newton's method from the parent's expression state
  (the default, much faster; see \[develop()\]).

## Value

An object of class \`grn_sim\`: a list with \`trajectory\` (data frame
of generation, trait means, phenotypic (co)variances \`P\_\*\`, mean
fitness and number of failed developments; plus \`M\_\*\` if
\`record_m\`, and \`gmean\_\*\`, \`G\_\*\` if \`record_genotype\`;
\`G3\_\*\` if \`record_third_moments\`), \`outcome\` (index of the peak
reached, or \`NA\`), \`generations\` run, and \`final\` (alleles,
genotypes, phenotypes and fitness of the last generation). With clonal
reproduction and no environmental noise the phenotypic covariance matrix
equals the genetic covariance matrix G.

## Details

Each locus has \`ploidy\` allele copies; the genotypic value of a locus
is the sum of its alleles. Mutations hit each allele copy with
probability \`mu\` per generation and add a normal deviate with sd
\`sigma\`; allelic values are kept within \`\[lower, upper\]\`.

With \`reproduction = "clonal"\` offspring copy their parent. With
\`"sexual"\`, each offspring has two parents (a female and a male if
\`sexes = TRUE\`), and inherits at each locus one randomly chosen allele
copy from each parent (free recombination; for haploids, the allele
comes from either parent).

The original analysis scripts correspond to \`ploidy = 1\`,
\`reproduction = "clonal"\`, \`selection = selection_relative_min()\`,
\`mu = 0.01\`, \`sigma = 0.05\`, \`N = 1000\`.

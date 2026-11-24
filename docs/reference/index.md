# Package index

## Networks

Genotype-phenotype maps built from gene regulatory networks.

- [`grn()`](https://diogro.github.io/grnland/reference/grn.md) :
  Hill-type gene regulatory network
- [`paper_motif()`](https://diogro.github.io/grnland/reference/paper_motif.md)
  : The six two-gene motifs of Netz, Melo & Milocco
- [`grn_custom()`](https://diogro.github.io/grnland/reference/grn_custom.md)
  : User-defined regulatory network
- [`gp_linear()`](https://diogro.github.io/grnland/reference/gp_linear.md)
  : Linear genotype-phenotype map
- [`loci_names()`](https://diogro.github.io/grnland/reference/loci_names.md)
  : Names of a model's loci

## Development and variation

From genotypes to phenotypes, and the variation mutation produces.

- [`default_genotype()`](https://diogro.github.io/grnland/reference/default_genotype.md)
  : Default genotype of a model (current values of the evolvable
  parameters)
- [`develop()`](https://diogro.github.io/grnland/reference/develop.md) :
  Develop genotypes into phenotypes
- [`phenotype()`](https://diogro.github.io/grnland/reference/phenotype.md)
  : Phenotype (trait values) of genotypes
- [`genotype_at()`](https://diogro.github.io/grnland/reference/genotype_at.md)
  : Genotype that develops into a given phenotype
- [`sensitivity()`](https://diogro.github.io/grnland/reference/sensitivity.md)
  : Derivatives of the developmental system at a steady state
- [`m_matrix()`](https://diogro.github.io/grnland/reference/m_matrix.md)
  : Mutational variance-covariance matrix (M-matrix)

## Landscapes

- [`landscape_perlin()`](https://diogro.github.io/grnland/reference/landscapes.md)
  [`landscape_gaussian()`](https://diogro.github.io/grnland/reference/landscapes.md)
  [`landscape_smoothmax()`](https://diogro.github.io/grnland/reference/landscapes.md)
  [`landscape_function()`](https://diogro.github.io/grnland/reference/landscapes.md)
  [`landscape_grid()`](https://diogro.github.io/grnland/reference/landscapes.md)
  : Adaptive landscapes
- [`find_peaks()`](https://diogro.github.io/grnland/reference/find_peaks.md)
  : Find local optima of a landscape
- [`lande_separatrices()`](https://diogro.github.io/grnland/reference/lande_separatrices.md)
  : Separatrices of the Lande-equation flow

## Simulation

- [`simulate_population()`](https://diogro.github.io/grnland/reference/simulate_population.md)
  : Simulate a population evolving on an adaptive landscape
- [`selection_relative_min()`](https://diogro.github.io/grnland/reference/selection.md)
  [`selection_proportional()`](https://diogro.github.io/grnland/reference/selection.md)
  [`selection_exponential()`](https://diogro.github.io/grnland/reference/selection.md)
  : Selection schemes
- [`stop_at_peaks()`](https://diogro.github.io/grnland/reference/stop_at_peaks.md)
  : Stop a simulation when the population reaches a peak
- [`nearest_peak()`](https://diogro.github.io/grnland/reference/nearest_peak.md)
  : Nearest peak to a phenotype
- [`develop_unique()`](https://diogro.github.io/grnland/reference/develop_unique.md)
  : Develop each distinct individual once
- [`dev_key_exact()`](https://diogro.github.io/grnland/reference/dev_key_exact.md)
  [`dev_key_genotype()`](https://diogro.github.io/grnland/reference/dev_key_exact.md)
  : Keys for de-duplicated development

## Package

- [`grnland`](https://diogro.github.io/grnland/reference/grnland-package.md)
  [`grnland-package`](https://diogro.github.io/grnland/reference/grnland-package.md)
  : grnland: Evolution of Gene Regulatory Networks on Adaptive
  Landscapes

# Selection schemes

Each returns a function mapping fitness values to parental sampling
weights. Individuals whose development failed (\`NA\` fitness) get
weight 0.

## Usage

``` r
selection_relative_min(eps = 0.001)

selection_proportional()

selection_exponential(s = 1)
```

## Arguments

- eps:

  Weight given to the least fit individual.

- s:

  Selection strength.

## Functions

- `selection_relative_min()`: Weights \\w - \min(w) + \epsilon\\. Used
  by the original analysis scripts. Selection is relative to the least
  fit individual, so it is very strong when the population has little
  fitness variance.

- `selection_proportional()`: Weights proportional to fitness
  (Wright-Fisher). Fitness must be non-negative.

- `selection_exponential()`: Weights \\\exp(s w)\\.

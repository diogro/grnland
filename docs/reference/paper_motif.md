# The six two-gene motifs of Netz, Melo & Milocco

Motifs from Gjuvsland et al. (2007) and Hether & Hohenlohe (2014), with
evolvable production rates \`alpha\`: 1: x1 -\> x2; 2: x1 -\| x2; 3: x1
-\> x2, x2 -\| x1; 4: x1 -\| x2, x2 -\| x1; 5: x1 -\> x2, x2 -\> x1; 6:
no interactions.

## Usage

``` r
paper_motif(k, theta = 5, gamma = 1, ...)
```

## Arguments

- k:

  Motif number, 1 to 6.

- theta, gamma:

  Threshold and decay rate (paper: 5 and 1).

- ...:

  Passed to \[grn()\], e.g. \`evolvable\`.

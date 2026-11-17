# Figure scripts

Scripts reproducing the figures of *Developmental Bias Shapes Evolution on
Complex Adaptive Landscapes* with the `grnland` package. They replace the
original scripts in the `GRN_Landscapes` repository.

| Script | Figures | Original |
|---|---|---|
| `fig1_linear_maps.R` | Fig 1 | none (simulation code was not in the repository) |
| `fig2_m_matrices.R` | Fig 2 | `Fig2.r` |
| `fig3_basins.R` | Fig 3, Fig S1, Fig S2 | `Figs34.r` |
| `fig4_fig5_complexity.R` | Fig 4, Fig 5 | `Fig5.r` (Fig 4); none for Fig 5 |
| `fig6_trajectories.R` | Fig 6 | `Fig6.r` |
| `figS3_perlin.R` | Fig S3 | none |
| `fig6_search.R` | (search for a Fig 6 example) | none |
| `motif_example.R` | `motif_example.png`: the six motifs reach different peaks from one start | none |

## Running

Run from this directory with BLAS threads capped:

```sh
cd analysis
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript fig3_basins.R
```

Environment variables:

- `QUICK=1` runs small grids and few replicates, writing to `output/quick` and
  `figures/quick`. Use it to check the pipeline.
- `N_CORES=k` sets the number of parallel workers (default: all cores but one).
- `REDO=1` ignores cached simulation results in `output/`.
- `SIGMA=0.05` (Figs 4 and 5 only) sets the mutational sd; the default is
  0.2, as in the manuscript.

Simulation results are cached as `.rds` files in `output/`; rerunning a script
only redraws the figures. All runs are seeded per job, so results do not
depend on the number of cores.

Approximate cost of the full runs (core-hours, at ~1 ms per generation):
Fig 3 ~200,000 runs and Fig 4/5 ~500,000 runs, each mostly a few hundred
generations; Fig 6 is 20 runs of N = 10,000 for 10,000 generations.

## Differences from the original scripts

Shared settings are in `common.R` (`PAPER`): haploid clonal populations,
N = 1000, mu = 0.01, sigma = 0.05, `w - min(w) + 0.001` selection, and
stopping when half the population is within ±0.1 of a peak.

- Phenotypes are always developed from the current genotype. This fixes the
  motif 5 bug in `Fig5.r`, which affected Figs 4 and 5.
- In Fig 4, starting points where any motif failed to converge are dropped for
  all motifs (the original dropped them only for motifs 1 and 6), and
  landscape seeds are reproducible.
- SVM 10-fold cross-validated accuracies are computed and saved
  (`output/fig3_svm_accuracy_*.csv`).
- Fig 2 uses the actual mutational variance, M = mu sigma^2 S S', with an
  explicit magnification, instead of Sigma = 0.1.
- Fig 6 runs 10 replicates per motif (originally 1) and reports how many reach
  each peak (`output/fig6_outcomes.csv`); M is not rescaled by an arbitrary
  factor, so panel C uses a log scale.
- Fig 1: the landscape seed, grid resolution and panel A starting points are
  not documented anywhere and are marked `TODO` in the script.

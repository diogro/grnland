# Grid screen of starting points, as in Fig 3 of the paper: for each motif,
# start monomorphic populations on a grid of phenotypes and record which peak
# each one reaches.
#
# Run with BLAS threads capped, one worker per core, e.g.
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript grid_screen.R

library(grnland)
library(parallel)

land <- landscape_grid(landscape_perlin(frequency = 0.5, seed = 1245), 0, 15)
peaks <- find_peaks(land, 0, 15, res = 500)

grid <- expand.grid(x1 = seq(2, 10, length.out = 41), x2 = seq(2, 10, length.out = 41))
jobs <- expand.grid(point = seq_len(nrow(grid)), motif = 1:6, rep = 1:10)

run_job <- function(k) {
  j <- jobs[k, ]
  start <- unlist(grid[j$point, ])
  sim <- simulate_population(paper_motif(j$motif), land, start = start,
                             N = 1000, generations = 10000,
                             stop = stop_at_peaks(peaks), record_every = 1000,
                             seed = k)
  last <- sim$trajectory[nrow(sim$trajectory), ]
  data.frame(motif = j$motif, rep = j$rep, x1 = start[1], x2 = start[2],
             peak = sim$outcome,
             nearest = nearest_peak(c(last$mean_x1, last$mean_x2), peaks),
             generations = sim$generations)
}

cores <- as.integer(Sys.getenv("N_CORES", detectCores()))
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores))
saveRDS(list(results = res, peaks = peaks), "grid_screen.rds")

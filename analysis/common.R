# Shared settings and helpers for the figure scripts.
# Run scripts from the analysis/ directory, e.g.
#   cd analysis && OPENBLAS_NUM_THREADS=1 Rscript fig3_basins.R
# Environment variables:
#   QUICK=1      small grids and few replicates, for testing the pipeline
#   N_CORES=k    number of parallel workers (default: all cores but one)
#   REDO=1       ignore cached simulation output and rerun

suppressPackageStartupMessages({
  library(grnland)
  library(tidyverse)
  library(parallel)
})

QUICK <- Sys.getenv("QUICK") == "1"
REDO <- Sys.getenv("REDO") == "1"
N_CORES <- as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 1)))
OUT_DIR <- if (QUICK) "output/quick" else "output"
FIG_DIR <- if (QUICK) "figures/quick" else "figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Simulation settings of the original analysis scripts. Note: the manuscript
# text says mu = 0.025 and diploid sexual populations; the code that produced
# the figures used these values.
PAPER <- list(N = 1000, generations = 10000, mu = 0.01, sigma = 0.05,
              ploidy = 1, reproduction = "clonal")

# Run fun(k) for k in seq_len(n) in parallel. Each job sets its own seed
# inside fun, so results do not depend on the number of cores.
run_jobs <- function(n, fun, cores = N_CORES) {
  res <- mclapply(seq_len(n), function(k) {
    tryCatch(fun(k), error = function(e) {
      message("job ", k, " failed: ", conditionMessage(e))
      NULL
    })
  }, mc.cores = cores, mc.preschedule = TRUE)
  bind_rows(res)
}

# Load cached output or compute and save it.
cached <- function(name, expr) {
  path <- file.path(OUT_DIR, paste0(name, ".rds"))
  if (file.exists(path) && !REDO) return(readRDS(path))
  t0 <- Sys.time()
  val <- expr
  message(sprintf("%s computed in %.1f min", name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(val, path)
  val
}

# One evolutionary run from a monomorphic start; returns the peak reached.
# Runs that do not meet the stopping rule within the generation limit are
# assigned the peak nearest to their final mean (as in the original scripts)
# and flagged with converged = FALSE.
run_to_peak <- function(model, land, start, peaks, seed, ...) {
  args <- modifyList(PAPER, list(...))
  sim <- do.call(simulate_population, c(list(model, land, start = start,
                                              stop = stop_at_peaks(peaks), seed = seed,
                                              record_every = args$generations), args))
  last <- sim$trajectory[nrow(sim$trajectory), ]
  end <- c(last$mean_x1, last$mean_x2)
  tibble(start_x1 = start[1], start_x2 = start[2],
         end_x1 = end[1], end_x2 = end[2],
         peak = if (is.na(sim$outcome)) nearest_peak(end, peaks) else sim$outcome,
         converged = !is.na(sim$outcome), generations = sim$generations)
}

# Basin borders from a radial SVM fitted to (start, outcome), plus the
# 10-fold cross-validated classification accuracy (%).
svm_borders <- function(d, res = 200, cost = 10, gamma = 0.5) {
  d <- d %>% transmute(x1 = start_x1, x2 = start_x2, peak = factor(peak))
  if (n_distinct(d$peak) < 2) return(list(borders = tibble(x1 = numeric(), x2 = numeric()), accuracy = 100))
  set.seed(1)  # reproducible cross-validation folds
  fit <- e1071::svm(peak ~ x1 + x2, data = d, kernel = "radial", cost = cost, gamma = gamma, cross = 10)
  xs <- seq(min(d$x1), max(d$x1), length.out = res)
  ys <- seq(min(d$x2), max(d$x2), length.out = res)
  g <- expand.grid(x1 = xs, x2 = ys)
  P <- matrix(as.integer(predict(fit, g)), res, res)
  B <- matrix(FALSE, res, res)
  B[-res, ] <- B[-res, ] | P[-1, ] != P[-res, ]
  B[, -res] <- B[, -res] | P[, -1] != P[, -res]
  idx <- which(B, arr.ind = TRUE)
  list(borders = tibble(x1 = xs[idx[, 1]], x2 = ys[idx[, 2]]), accuracy = fit$tot.accuracy)
}

# Landscape values on a grid, for plotting.
landscape_df <- function(land, lower = 0, upper = 11.5, res = 300) {
  g <- expand.grid(x = seq(lower, upper, length.out = res), y = seq(lower, upper, length.out = res))
  g$z <- land(as.matrix(g))
  g
}

# Covariance ellipse (level-quantile contour of a bivariate normal).
ellipse_path <- function(M, centre, level = 0.95, n = 60) {
  e <- eigen(M, symmetric = TRUE)
  if (any(e$values <= 0)) return(NULL)
  t <- seq(0, 2 * pi, length.out = n)
  r <- sqrt(qchisq(level, 2))
  pts <- r * cbind(cos(t), sin(t)) %*% diag(sqrt(e$values)) %*% t(e$vectors)
  tibble(ex = pts[, 1] + centre[1], ey = pts[, 2] + centre[2])
}

motif_label <- function(x) paste("motif", x)

theme_set(theme_minimal())

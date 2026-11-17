# Fig 1: evolution under three linear genotype-phenotype maps x = B alpha.
#   A: single Gaussian peak, trajectories from four starting points.
#   B-D: two-peak Perlin landscape, outcomes from a grid of starting points,
#        basin border inferred by SVM, plus the Lande-equation separatrix
#        (G proportional to B B') for comparison.
#
# The simulation code for this figure is not in the original repository
# (GRN_Landscapes/Fig1.py draws a different, analytical figure). Settings not
# stated in the manuscript are marked TODO: replace them with the original
# values if they can be recovered.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

B <- list(B1 = diag(2),
          B2 = matrix(c(5, 4, 4, 5) / 3, 2),
          B3 = matrix(c(5, -4, -4, 5) / 3, 2))
models <- lapply(B, gp_linear)
# Allelic values may be negative for linear maps (B2 and B3 need negative
# alpha in parts of the region).
LOWER <- -Inf

# ---- A: single peak ----------------------------------------------------------
single <- landscape_gaussian(c(5.5, 5.5), omega = 3^2)   # from the comment in results.tex
starts_A <- rbind(c(1, 4), c(4, 10), c(10, 7), c(7, 1))  # TODO: original starting points unknown
trajA <- cached("fig1_single", run_jobs(nrow(starts_A) * 3, function(k) {
  b <- (k - 1) %/% nrow(starts_A) + 1
  s <- (k - 1) %% nrow(starts_A) + 1
  sim <- simulate_population(models[[b]], single, start = starts_A[s, ], N = PAPER$N,
                             generations = PAPER$generations, mu = PAPER$mu, sigma = PAPER$sigma,
                             lower = LOWER, record_every = 10, seed = 10 + k,
                             stop = stop_at_peaks(matrix(c(5.5, 5.5), 1)))
  as_tibble(sim$trajectory) %>% mutate(map = names(B)[b], start = s)
}))

# ---- B-D: two peaks ----------------------------------------------------------
# TODO: original landscape unknown. Pick the first seed whose frequency-0.1
# Perlin landscape has exactly two optima in [0, 11]^2, both inside [1, 9]^2.
FREQ <- 0.1
seed <- 1
repeat {
  perlin <- landscape_perlin(frequency = FREQ, seed = seed)
  pk <- find_peaks(perlin, 0, 15, res = 300)
  near <- pk$x1 < 11 & pk$x2 < 11
  inner <- pk$x1 > 1 & pk$x1 < 9 & pk$x2 > 1 & pk$x2 < 9
  if (sum(near) == 2 && sum(inner) == 2) break
  seed <- seed + 1
}
message("Fig 1 two-peak landscape: frequency ", FREQ, ", seed ", seed)
land <- landscape_grid(perlin, 0, 15)
peaks <- find_peaks(perlin, 0, 15, res = 500)

GRID_N <- if (QUICK) 7 else 25   # TODO: original grid resolution unknown
REPS <- if (QUICK) 2 else 10
grid <- expand.grid(x1 = seq(2, 8, length.out = GRID_N), x2 = seq(0, 6, length.out = GRID_N))
jobs <- expand.grid(point = seq_len(nrow(grid)), map = seq_along(B), rep = seq_len(REPS))
res <- cached("fig1_two_peaks", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  run_to_peak(models[[j$map]], land, unlist(grid[j$point, ]), peaks,
              seed = 1000 + k, lower = LOWER) %>%
    mutate(map = names(B)[j$map], rep = j$rep)
}))

svm <- lapply(names(B), function(b) svm_borders(filter(res, map == b)))
borders <- bind_rows(lapply(seq_along(B), function(i) mutate(svm[[i]]$borders, map = names(B)[i])))
print(tibble(map = names(B), cv_accuracy = sapply(svm, `[[`, "accuracy")))
lande <- bind_rows(lapply(names(B), function(b)
  mutate(lande_separatrices(models[[b]], perlin, lower = -5, upper = 15, step = 0.02), map = b)))

M_text <- sapply(names(B), function(b) {
  M <- m_matrix(models[[b]], c(1, 1), mu = PAPER$mu, sigma = PAPER$sigma)
  sprintf("M x 1e5 = [%.1f, %.1f; %.1f, %.1f]", 1e5 * M[1, 1], 1e5 * M[1, 2], 1e5 * M[2, 1], 1e5 * M[2, 2])
})

bgA <- landscape_df(single, 0, 11, res = 200)
pA <- ggplot(bgA, aes(x, y)) + geom_raster(aes(fill = z)) + scale_fill_viridis_c(name = "fitness") +
  geom_path(data = trajA, aes(mean_x1, mean_x2, colour = map, group = interaction(map, start)), linewidth = 0.6) +
  scale_colour_manual(values = c(B1 = "red", B2 = "white", B3 = "orange")) +
  theme(legend.key = element_rect(fill = "grey40", colour = NA)) +
  coord_equal() + labs(x = expression(x[1]), y = expression(x[2]), colour = "map")

agree <- res %>% count(map, start_x1, start_x2, peak) %>% group_by(map, start_x1, start_x2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(freq = n / REPS)
ends <- res %>% mutate(px = peaks$x1[peak], py = peaks$x2[peak])
bg <- landscape_df(land, 0, 10)
pBD <- ggplot(bg, aes(x, y)) + geom_raster(aes(fill = z)) + scale_fill_viridis_c(name = "fitness") +
  geom_segment(data = ends, aes(x = start_x1, y = start_x2, xend = px, yend = py), inherit.aes = FALSE,
               colour = "orange", linewidth = 0.1) +
  ggnewscale::new_scale_fill() +
  geom_point(data = agree, aes(start_x1, start_x2, colour = freq), size = 0.4, inherit.aes = FALSE) +
  scale_colour_gradient(low = "white", high = "black", name = "max peak frequency") +
  geom_tile(data = borders, aes(x1, x2), fill = "white", inherit.aes = FALSE) +
  geom_path(data = lande, aes(x1, x2, group = branch), colour = "black", linetype = 2,
            linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(data = filter(peaks, x1 < 10, x2 < 10), aes(x1, x2), colour = "red", shape = 4,
             size = 2, inherit.aes = FALSE) +
  coord_equal(xlim = c(0, 10), ylim = c(0, 10)) +
  facet_wrap(~map, labeller = labeller(map = function(b) paste(b, M_text[b], sep = "\n"))) +
  labs(x = expression(x[1]), y = expression(x[2]))

ggsave(file.path(FIG_DIR, "fig1.png"),
       (pA + labs(tag = "A")) + (pBD + labs(tag = "B")) + plot_layout(widths = c(1, 3)),
       width = 15, height = 5.5, bg = "white")

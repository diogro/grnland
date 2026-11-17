# Fig 6: evolutionary trajectories of motifs 1 and 5 from (nearly) the same
# M-matrix, on a landscape with two adjacent peaks. Tracks mean phenotype, the
# G-matrix (phenotypic covariance; equal to G for clonal populations without
# environmental noise) and the M-matrix at the mean genotype.
# Port of GRN_Landscapes/Fig6.r.
#
# Changes from the original:
#  - REPS replicates per motif (original: 1); panels A-C show replicate 1 and
#    the proportion of replicates reaching each peak is reported;
#  - M is computed with the actual mutational variance (mu * sigma^2); the
#    original used Sigma = 0.1 and rescaled M by an arbitrary 0.001 for panel C.
#    Panel C therefore shows variances on a log scale, and panel B draws M
#    rescaled to the trace of G so that shapes can be compared;
#  - the steepest-ascent path is computed directly from the landscape.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

START <- c(1.12, 6.4)
REPS <- if (QUICK) 2 else 10
N <- if (QUICK) 1000 else 10000
GENS <- if (QUICK) 2000 else 10000

# Two shifted, squared Perlin layers combined by a smooth maximum.
layer <- function(offset) landscape_perlin(frequency = 0.05, seed = 76, offset = offset,
                                           transform = function(z) ((z + 1) / 2)^2)
perlin <- landscape_smoothmax(layer(c(0, 1)), layer(c(-1, 0)), beta = 100)
land <- landscape_grid(perlin, 0, 15)
peaks <- find_peaks(perlin, 0, 15, res = 500)

jobs <- expand.grid(rep = seq_len(REPS), motif = c(1, 5))
res <- cached("fig6", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  sim <- simulate_population(paper_motif(j$motif), land, start = START,
                             N = N, generations = GENS, mu = PAPER$mu, sigma = PAPER$sigma,
                             record_every = 10, record_m = TRUE, seed = 600 + k)
  last <- sim$trajectory[nrow(sim$trajectory), ]
  sim$trajectory %>% as_tibble() %>%
    mutate(motif = j$motif, rep = j$rep,
           final_peak = nearest_peak(c(last$mean_x1, last$mean_x2), peaks))
}, cores = min(N_CORES, nrow(jobs))))

outcomes <- res %>% distinct(motif, rep, final_peak) %>% count(motif, final_peak)
print(outcomes)
write_csv(outcomes, file.path(OUT_DIR, "fig6_outcomes.csv"))

tr <- filter(res, rep == 1)

# ---- A: landscape, trajectories, steepest ascent ----------------------------
ascent <- function(z, ds = 0.01, n = 3000, h = 1e-4) {
  path <- matrix(NA, n, 2)
  for (i in seq_len(n)) {
    path[i, ] <- z
    g <- (perlin(rbind(z + c(h, 0), z + c(0, h))) - perlin(rbind(z - c(h, 0), z - c(0, h)))) / (2 * h)
    if (sqrt(sum(g^2)) < 1e-6) break
    z <- z + ds * g / sqrt(sum(g^2))
  }
  as_tibble(path[complete.cases(path), ], .name_repair = ~ c("x", "y"))
}
asc <- ascent(START)
bg <- landscape_df(land, 0, 15)
pA <- ggplot(bg, aes(x, y)) +
  geom_raster(aes(fill = z)) + scale_fill_viridis_c(name = "fitness") +
  geom_path(data = asc, colour = "red", linewidth = 1) +
  geom_point(data = tr, aes(mean_x1, mean_x2, colour = generation), size = 0.5) +
  scale_colour_gradient(low = "darkblue", high = "lightblue") +
  geom_point(data = peaks, aes(x1, x2), shape = 4, size = 1.5, stroke = 1, inherit.aes = FALSE) +
  coord_equal() + facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
  labs(x = expression(x[1]), y = expression(x[2]))

# ---- B: trajectories with G and M ellipses ----------------------------------
mats <- tr %>% filter(generation %% 200 == 0) %>%
  mutate(G = pmap(list(P_x1_x1, P_x1_x2, P_x2_x2), ~ matrix(c(..1, ..2, ..2, ..3), 2)),
         M = pmap(list(M_x1_x1, M_x1_x2, M_x2_x2), ~ matrix(c(..1, ..2, ..2, ..3), 2)),
         M = map2(M, G, ~ .x * sum(diag(.y)) / sum(diag(.x))))
ells <- function(col) mats %>%
  mutate(e = map2(.data[[col]], map2(mean_x1, mean_x2, c), ~ ellipse_path(2^2 * .x, .y))) %>%
  select(motif, generation, e) %>% unnest(e)
pB <- ggplot(tr, aes(mean_x1, mean_x2)) +
  geom_path() +
  geom_path(data = ells("M"), aes(ex, ey, group = generation), colour = "#00BA38", alpha = 0.6) +
  geom_path(data = ells("G"), aes(ex, ey, group = generation), colour = "#F8766D", alpha = 0.6) +
  scale_x_continuous(n.breaks = 3) +
  coord_equal() + facet_wrap(~motif) +
  labs(x = expression(x[1]), y = expression(x[2])) + theme(strip.text = element_blank())

# ---- C: eccentricity, angle, total variance ---------------------------------
summ <- tr %>%
  pivot_longer(c(P_x1_x1, P_x1_x2, P_x2_x2, M_x1_x1, M_x1_x2, M_x2_x2),
               names_to = c("matrix", ".value"), names_pattern = "([PM])_(.*)") %>%
  mutate(matrix = if_else(matrix == "P", "G-matrix", "M-matrix"),
         disc = sqrt((x1_x1 - x2_x2)^2 + 4 * x1_x2^2),
         lmax = (x1_x1 + x2_x2 + disc) / 2, lmin = (x1_x1 + x2_x2 - disc) / 2,
         eccentricity = sqrt(1 - pmax(lmin, 0) / lmax),
         angle = 0.5 * atan2(2 * x1_x2, x1_x1 - x2_x2) * 180 / pi,
         `total variance` = x1_x1 + x2_x2) %>%
  filter(generation > 0)
panel <- function(var, log = FALSE) {
  p <- ggplot(summ, aes(generation, .data[[var]], colour = matrix)) + geom_line() +
    facet_wrap(~motif, labeller = labeller(motif = motif_label)) + labs(y = var, colour = NULL)
  if (log) p + scale_y_log10() else p
}
pC <- (panel("eccentricity") / panel("angle") / panel("total variance", log = TRUE)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig6.png"),
       ((pA + labs(tag = "A")) / (pB + labs(tag = "B")) | (pC + plot_annotation(tag_levels = list("C")))),
       width = 11, height = 7, bg = "white")

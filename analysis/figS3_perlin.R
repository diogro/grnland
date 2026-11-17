# Fig S3: Perlin landscapes of increasing frequency (same seed) and their
# spatial autocorrelation. Not in the original repository; implemented from
# the manuscript's description.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

FREQS <- c(0.1, 0.3, 0.6, 1.2, 1.8)
SEED <- 1234

maps <- bind_rows(lapply(FREQS, function(f) {
  landscape_df(landscape_perlin(frequency = f, seed = SEED, offset = 0.5), 0, 10, res = 250) %>%
    mutate(frequency = f)
}))
pA <- ggplot(maps, aes(x, y, fill = z)) + geom_raster() +
  scale_fill_viridis_c(name = "fitness") + coord_equal() +
  facet_wrap(~frequency, nrow = 1, labeller = labeller(frequency = function(f) paste("frequency", f))) +
  labs(x = expression(x[1]), y = expression(x[2]))

# Correlation of fitness between points a distance d apart, in random
# directions, over a large region.
set.seed(1)
n <- 20000
dists <- seq(0, 10, by = 0.1)
acf <- bind_rows(lapply(FREQS, function(f) {
  land <- landscape_perlin(frequency = f, seed = SEED, offset = 0.5)
  X <- matrix(runif(2 * n, 0, 50), ncol = 2)
  w0 <- land(X)
  tibble(frequency = f, distance = dists, correlation = sapply(dists, function(d) {
    a <- runif(n, 0, 2 * pi)
    cor(w0, land(X + d * cbind(cos(a), sin(a))))
  }))
}))
pB <- ggplot(acf, aes(distance, correlation, colour = factor(frequency))) + geom_line() +
  geom_hline(yintercept = 0, linetype = 2) + labs(colour = "frequency")

ggsave(file.path(FIG_DIR, "figS3.png"), pA / pB + plot_annotation(tag_levels = "A"),
       width = 10, height = 5.5, bg = "white")

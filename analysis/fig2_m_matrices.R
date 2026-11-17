# Fig 2: M-matrices and sensitivity vectors of the six motifs over phenotype
# space. Port of GRN_Landscapes/Fig2.r.
#
# Change from the original: the original used Sigma = diag(2 * sigma) = 0.1,
# which is not the mutational variance. Here M = mu * sigma^2 * S S' (haploid,
# the simulation settings) and the ellipses are magnified by an explicit
# factor, reported in the output.

source("common.R")

MAG <- 40   # magnification of ellipse axes (standard deviations)
pts <- expand.grid(x1 = seq(2, 10, length.out = 5), x2 = seq(2, 10, length.out = 5))

rows <- list()
for (k in 1:6) {
  m <- paper_motif(k)
  for (r in seq_len(nrow(pts))) {
    x <- unlist(pts[r, ])
    g <- genotype_at(m, x)
    S <- sensitivity(m, g)$S
    M <- m_matrix(m, g, mu = PAPER$mu, sigma = PAPER$sigma, ploidy = PAPER$ploidy)
    rows[[length(rows) + 1]] <- tibble(
      motif = k, x1 = x[1], x2 = x[2],
      s1_x = S[1, 1], s1_y = S[2, 1], s2_x = S[1, 2], s2_y = S[2, 2],
      ell = list(ellipse_path(MAG^2 * M, x)))
  }
}
d <- bind_rows(rows)
ell <- d %>% select(motif, x1, x2, ell) %>% unnest(ell)

p <- ggplot(ell, aes(ex, ey, group = interaction(motif, x1, x2))) +
  geom_path() +
  geom_segment(data = d, aes(x = x1, y = x2, xend = x1 + s1_x, yend = x2 + s1_y), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour = "blue", linewidth = 0.4) +
  geom_segment(data = d, aes(x = x1, y = x2, xend = x1 + s2_x, yend = x2 + s2_y), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour = "red", linewidth = 0.4) +
  coord_equal() +
  facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(strip.text = element_text(face = "bold"), panel.spacing.x = unit(3, "lines"))

ggsave(file.path(FIG_DIR, "fig2.png"), p, width = 7, height = 4, bg = "white")
message(sprintf("Fig 2: 95%% ellipses of M (mu = %g, sigma = %g, ploidy %d), axes magnified %gx",
                PAPER$mu, PAPER$sigma, PAPER$ploidy, MAG))

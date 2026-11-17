# Fig 4 and Fig 5: effect of landscape complexity on evolutionary outcomes.
# For each Perlin frequency, 10 random landscapes; every motif is run once from
# each point of a grid of starting phenotypes.
#   Fig 4A: outcome mismatch with motif 6 (no developmental bias) vs frequency.
#   Fig 4B: proportion of identical outcomes between all pairs of motifs.
#   Fig 5:  angle between the direction of evolutionary divergence (start ->
#           peak reached) and the leading eigenvector of M at the end point.
# Port of GRN_Landscapes/Fig5.r (Fig 4); Fig 5 had no code in the original
# repository and is implemented here from the manuscript's description.
#
# Changes from the original:
#  - motif 5 phenotypes are always developed from the mutated genotype (the
#    original kept the pre-mutation phenotype; see the package README);
#  - landscape seeds are drawn reproducibly;
#  - starting points where any motif failed to converge are dropped for all
#    comparisons (the original dropped them only for motifs 1 and 6).
#
# SIGMA (env var) sets the mutational sd; the manuscript uses 0.2 here and
# mentions a 0.05 run.

source("common.R")

SIGMA <- as.numeric(Sys.getenv("SIGMA", "0.2"))
FREQS <- if (QUICK) c(0.1, 0.6, 1.8) else c(0.1, 0.3, 0.6, 1.2, 1.8)
N_LAND <- if (QUICK) 2 else 10
GRID_N <- if (QUICK) 9 else 41
set.seed(20261117)
lands <- crossing(frequency = FREQS, landscape = seq_len(N_LAND)) %>%
  mutate(seed = sample.int(10000, n()))
grid <- expand.grid(x1 = seq(2, 10, length.out = GRID_N), x2 = seq(2, 10, length.out = GRID_N))
tag <- paste0("sigma", SIGMA)

make_land <- function(l) {
  # offset by 0.5 to keep the grid of starting points off the Perlin lattice
  perlin <- landscape_perlin(frequency = l$frequency, seed = l$seed, offset = 0.5)
  list(perlin = perlin, land = landscape_grid(perlin, 0, 30, step = 0.01),
       peaks = find_peaks(perlin, 0, 30, res = 1000))
}

res <- cached(paste0("fig4_", tag), {
  bind_rows(lapply(seq_len(nrow(lands)), function(L) {
    l <- lands[L, ]
    message(sprintf("frequency %g, landscape %d", l$frequency, l$landscape))
    ld <- make_land(l)
    jobs <- expand.grid(point = seq_len(nrow(grid)), motif = 1:6)
    run_jobs(nrow(jobs), function(k) {
      j <- jobs[k, ]
      run_to_peak(paper_motif(j$motif), ld$land, unlist(grid[j$point, ]), ld$peaks,
                  seed = 1e6 * L + k, sigma = SIGMA) %>%
        mutate(motif = j$motif)
    }) %>%
      mutate(frequency = l$frequency, landscape = l$landscape,
             peak_x1 = ld$peaks$x1[peak], peak_x2 = ld$peaks$x2[peak])
  }))
})
message(sprintf("%d runs, %d did not converge", nrow(res), sum(!res$converged)))

# ---- Fig 4: outcome mismatch -------------------------------------------------
wide <- res %>%
  group_by(frequency, landscape, start_x1, start_x2) %>%
  filter(all(converged)) %>%
  ungroup() %>%
  select(frequency, landscape, start_x1, start_x2, motif, peak) %>%
  pivot_wider(names_from = motif, values_from = peak, names_prefix = "m")

pairs <- combn(1:6, 2, simplify = FALSE)
pairwise <- bind_rows(lapply(pairs, function(p) {
  wide %>% group_by(frequency, landscape) %>%
    summarise(same = mean(.data[[paste0("m", p[1])]] == .data[[paste0("m", p[2])]]), .groups = "drop") %>%
    mutate(motif_1 = p[1], motif_2 = p[2])
}))
write_csv(pairwise, file.path(OUT_DIR, paste0("fig4_pairwise_", tag, ".csv")))

pA <- pairwise %>% filter(motif_2 == 6) %>%
  ggplot(aes(frequency, 1 - same)) +
  geom_point(alpha = 0.5, size = 1) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(breaks = FREQS) +
  facet_wrap(~motif_1, labeller = labeller(motif_1 = motif_label), nrow = 1) +
  labs(x = "landscape complexity (Perlin frequency)", y = "outcome mismatch with motif 6")

pooled <- pairwise %>% group_by(motif_1, motif_2) %>% summarise(same = mean(same), .groups = "drop")
pB <- if (requireNamespace("ggraph", quietly = TRUE)) {
  g <- igraph::graph_from_data_frame(
    pooled %>% transmute(from = motif_label(motif_1), to = motif_label(motif_2), same),
    directed = FALSE)
  set.seed(1)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = same, alpha = same)) +
    ggraph::geom_node_point(size = 5) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE) +
    ggraph::scale_edge_width(range = c(0.5, 3), name = "same outcomes") +
    ggraph::scale_edge_alpha(range = c(0.3, 1), name = "same outcomes") +
    theme_void()
} else {
  ggplot(pooled, aes(factor(motif_1), factor(motif_2), fill = same)) + geom_tile() +
    labs(x = "motif", y = "motif", fill = "same outcomes")
}

ggsave(file.path(FIG_DIR, paste0("fig4_", tag, ".png")),
       patchwork::wrap_plots(pA + labs(tag = "A"), pB + labs(tag = "B"), ncol = 1, heights = c(1, 1.4)),
       width = 9, height = 8, bg = "white")

# ---- Fig 5: alignment of divergence with M ----------------------------------
# M at each peak reached, per motif (cached: many runs share peaks).
peak_M <- res %>% filter(converged) %>% distinct(motif, peak_x1, peak_x2)
peak_M$angle_vec <- mapply(function(k, x1, x2) {
  M <- tryCatch(m_matrix(paper_motif(k), phenotype = c(x1, x2)), error = function(e) NULL)
  if (is.null(M)) return(NA_real_)
  v <- eigen(M, symmetric = TRUE)$vectors[, 1]
  atan2(v[2], v[1])
}, peak_M$motif, peak_M$peak_x1, peak_M$peak_x2)

align <- res %>% filter(converged) %>%
  left_join(peak_M, by = c("motif", "peak_x1", "peak_x2")) %>%
  mutate(dx = peak_x1 - start_x1, dy = peak_x2 - start_x2) %>%
  filter(sqrt(dx^2 + dy^2) > 0.1, !is.na(angle_vec)) %>%
  mutate(angle = acos(abs(cos(atan2(dy, dx) - angle_vec))) * 180 / pi)
# Note: for motif 6, M is exactly isotropic, so its "leading eigenvector" is
# an arbitrary axis returned by eigen(); the angles then only reflect the
# distribution of divergence directions.

p5 <- ggplot(align, aes(factor(frequency), angle)) +
  geom_violin(fill = "grey85", colour = NA) +
  stat_summary(fun = median, geom = "point") +
  geom_hline(yintercept = 45, linetype = 2) +
  facet_wrap(~motif, labeller = labeller(motif = motif_label), nrow = 1) +
  labs(x = "landscape complexity (Perlin frequency)",
       y = expression("angle between divergence and " * M[max] * " (degrees)"))
ggsave(file.path(FIG_DIR, paste0("fig5_", tag, ".png")), p5, width = 10, height = 3.5, bg = "white")

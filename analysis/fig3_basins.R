# Fig 3, Fig S1, Fig S2: evolutionary outcomes from a grid of starting points
# for all six motifs on a simple (Perlin frequency 0.1) and a complex
# (frequency 0.5) landscape, basin borders inferred by SVM, and separatrices
# predicted by the Lande equation with G proportional to M.
# Port of GRN_Landscapes/Figs34.r.
#
# Changes from the original: SVM cross-validated accuracy is computed (it was
# a to-do); non-converged runs are flagged and counted.

source("common.R")
suppressPackageStartupMessages(library(ggnewscale))

LANDSCAPES <- tibble(name = c("simple", "complex"), frequency = c(0.1, 0.5), seed = c(1234, 1245))
GRID_N <- if (QUICK) 9 else 41
REPS <- if (QUICK) 2 else 10
grid <- expand.grid(x1 = seq(2, 10, length.out = GRID_N), x2 = seq(2, 10, length.out = GRID_N))

for (L in seq_len(nrow(LANDSCAPES))) {
  ls <- LANDSCAPES[L, ]
  perlin <- landscape_perlin(frequency = ls$frequency, seed = ls$seed)
  land <- landscape_grid(perlin, 0, 15)
  peaks <- find_peaks(perlin, 0, 15, res = 500)

  jobs <- expand.grid(point = seq_len(nrow(grid)), motif = 1:6, rep = seq_len(REPS))
  res <- cached(paste0("fig3_", ls$name), run_jobs(nrow(jobs), function(k) {
    j <- jobs[k, ]
    run_to_peak(paper_motif(j$motif), land, unlist(grid[j$point, ]), peaks,
                seed = 100000 * L + k) %>%
      mutate(motif = j$motif, rep = j$rep)
  }))
  message(sprintf("%s landscape: %d runs, %d did not converge", ls$name, nrow(res), sum(!res$converged)))

  svm <- lapply(1:6, function(k) svm_borders(filter(res, motif == k)))
  borders <- bind_rows(lapply(1:6, function(k) mutate(svm[[k]]$borders, motif = k)))
  accuracy <- tibble(landscape = ls$name, motif = 1:6, cv_accuracy = sapply(svm, `[[`, "accuracy"))
  print(accuracy)
  write_csv(accuracy, file.path(OUT_DIR, paste0("fig3_svm_accuracy_", ls$name, ".csv")))

  lande <- cached(paste0("fig3_lande_", ls$name), bind_rows(mclapply(1:6, function(k) {
    mutate(lande_separatrices(paper_motif(k), perlin, lower = 0, upper = 12, step = 0.01), motif = k)
  }, mc.cores = min(6, N_CORES)))) %>% filter(between(x1, 2, 10), between(x2, 2, 10))

  # frequency of the most common outcome at each starting point
  agree <- res %>% count(motif, start_x1, start_x2, peak) %>%
    group_by(motif, start_x1, start_x2) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% mutate(freq = n / REPS)
  ends <- res %>% mutate(px = peaks$x1[peak], py = peaks$x2[peak])
  bg <- landscape_df(land, 0, 11.5)

  outcome_plot <- function(motifs, ref_border = TRUE, lande_lines = FALSE) {
    p <- ggplot(bg, aes(x, y)) +
      geom_raster(aes(fill = z)) +
      viridis::scale_fill_viridis(name = "fitness") +
      geom_segment(data = filter(ends, motif %in% motifs),
                   aes(x = start_x1, y = start_x2, xend = px, yend = py), inherit.aes = FALSE,
                   arrow = arrow(length = unit(0.05, "cm")), colour = "orange", linewidth = 0.1) +
      geom_point(data = filter(ends, motif %in% motifs), aes(px, py), inherit.aes = FALSE,
                 shape = 4, size = 1, colour = "red")
    if (ref_border) {
      ref <- filter(borders, motif == 6) %>% select(-motif) %>% crossing(motif = motifs)
      p <- p + geom_tile(data = ref, aes(x1, x2), fill = "red", inherit.aes = FALSE)
    }
    p <- p + new_scale_fill() +
      geom_point(data = filter(agree, motif %in% motifs), aes(start_x1, start_x2, colour = freq),
                 size = 0.2, inherit.aes = FALSE) +
      scale_colour_gradient(low = "white", high = "black", name = "max peak frequency") +
      geom_tile(data = filter(borders, motif %in% motifs), aes(x1, x2), fill = "white", inherit.aes = FALSE)
    if (lande_lines)
      p <- p + geom_path(data = filter(lande, motif %in% motifs), aes(x1, x2, group = branch),
                         colour = "black", linewidth = 0.2, inherit.aes = FALSE)
    p + coord_equal(xlim = c(0, 11.5), ylim = c(0, 11.5)) +
      facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
      labs(x = expression(x[1]), y = expression(x[2]))
  }

  ggsave(file.path(FIG_DIR, paste0("fig3_", ls$name, ".png")), outcome_plot(c(1, 2)),
         width = 8, height = 4, bg = "white")
  ggsave(file.path(FIG_DIR, paste0("figS", L, ".png")), outcome_plot(1:6, ref_border = FALSE, lande_lines = TRUE),
         width = 10, height = 7, bg = "white")
}

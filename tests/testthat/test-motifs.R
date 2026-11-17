grid <- expand.grid(a1 = c(0.5, 2, 6, 11), a2 = c(0.7, 3, 8, 14))

test_that("numerical steady states match the closed forms of the paper motifs", {
  for (k in 1:6) {
    m <- paper_motif(k)
    G <- as.matrix(grid)
    expected <- t(apply(G, 1, function(g) legacy_GtoP[[k]](g[1], g[2])))
    # motif 5: the closed form is the positive state, which exists only if a1*a2 > theta^2
    keep <- if (k == 5) G[, 1] * G[, 2] > 25 else rep(TRUE, nrow(G))
    d <- develop(m, G[keep, , drop = FALSE], x_start = c(1, 1))
    expect_true(all(d$status == 0), info = paste("motif", k))
    expect_equal(unname(d$X), unname(expected[keep, , drop = FALSE]), tolerance = 1e-8,
                 info = paste("motif", k))
  }
})

test_that("motif 5 goes to the zero state when no positive state exists", {
  d <- develop(paper_motif(5), c(2, 3), x_start = c(1, 1))
  expect_equal(d$status, 0L)
  expect_equal(unname(d$X[1, ]), c(0, 0), tolerance = 1e-8)
})

test_that("analytical sensitivities match the closed forms of Fig2.r", {
  pts <- expand.grid(x1 = c(2, 4.5, 8), x2 = c(2, 6, 10))
  for (k in 1:6) {
    m <- paper_motif(k)
    for (r in seq_len(nrow(pts))) {
      x <- unlist(pts[r, ])
      g <- genotype_at(m, x)
      S <- sensitivity(m, g)$S
      expect_equal(unname(S), legacy_S[[k]](x[1], x[2]), tolerance = 1e-7,
                   info = sprintf("motif %d at (%g, %g)", k, x[1], x[2]))
    }
  }
})

test_that("genotype_at inverts development", {
  for (k in 1:6) {
    m <- paper_motif(k)
    x <- c(3.3, 7.1)
    g <- genotype_at(m, x)
    expect_equal(unname(phenotype(m, g, x_start = x)), x, tolerance = 1e-8)
    expect_equal(g, unname(motifPG(k, x)), tolerance = 1e-10)
  }
})

test_that("numerical inversion works when loci are not production rates", {
  m <- grn(c("x1 -> x2", "x2 -| x1"), alpha = c(4, 6), theta = 5,
           evolvable = c("theta"))
  g <- genotype_at(m, c(2, 3), start = c(5, 5))
  expect_equal(unname(phenotype(m, g, x_start = c(2, 3))), c(2, 3), tolerance = 1e-8)
})

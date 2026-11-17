test_that("edge notation and matrix notation give the same model", {
  a <- grn(c("x1 -> x2", "x2 -| x1"))
  reg <- matrix(c(0, 1, -1, 0), 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
  b <- grn(reg)
  expect_identical(a$regulation, b$regulation)
  expect_error(grn("x1 => x2"), "cannot parse")
})

test_that("evolvable parameter specs expand correctly", {
  m <- grn(c("a -> b", "b -| c", "c -> a"), evolvable = c("alpha", "theta", "hill[b,a]"))
  expect_equal(loci_names(m), c("alpha[a]", "alpha[b]", "alpha[c]",
                                "theta[a,c]", "theta[b,a]", "theta[c,b]", "hill[b,a]"))
})

test_that("analytical derivatives agree with finite differences for all parameter types", {
  m <- grn(c("a -> b", "b -| c", "c -> a", "a -| c"), alpha = c(3, 4, 5), gamma = c(1, 0.8, 1.2),
           theta = 2, hill = 2, evolvable = c("alpha", "gamma", "theta", "hill"))
  g <- default_genotype(m)
  s <- sensitivity(m, g)
  num <- sapply(seq_along(g), function(l) {
    h <- 1e-5 * max(1, g[l])
    gp <- g; gm <- g
    gp[l] <- gp[l] + h; gm[l] <- gm[l] - h
    (phenotype(m, gp, x_start = s$x) - phenotype(m, gm, x_start = s$x)) / (2 * h)
  })
  expect_equal(unname(s$S), unname(num), tolerance = 1e-5)
})

test_that("custom models reproduce the Hill model", {
  hill <- paper_motif(3)
  custom <- grn_custom(function(x, p) c(p[["a1"]] * (1 - x[2] / (5 + x[2])) - x[1],
                                        p[["a2"]] * x[1] / (5 + x[1]) - x[2]),
                       params = c(a1 = 1, a2 = 1), genes = c("x1", "x2"))
  g <- c(6, 9)
  expect_equal(phenotype(custom, g), phenotype(hill, g), tolerance = 1e-8)
  expect_equal(unname(m_matrix(custom, g)), unname(m_matrix(hill, g)), tolerance = 1e-6)
})

test_that("M-matrix of a linear map is ploidy * mu * sigma^2 * B B'", {
  B <- matrix(c(5, 4, 4, 5) / 3, 2)
  lin <- grn_custom(function(x, p) as.numeric(B %*% c(p[["a1"]], p[["a2"]])) - x,
                    params = c(a1 = 3, a2 = 3), genes = 2)
  M <- m_matrix(lin, c(3, 3), mu = 0.025, sigma = 0.05, ploidy = 2)
  expect_equal(unname(M), 2 * 0.025 * 0.05^2 * B %*% t(B), tolerance = 1e-8)
})

test_that("trait subsets restrict the M-matrix and phenotype", {
  m <- grn(c("a -> b", "b -> c"), theta = 2, alpha = 3, traits = c("a", "c"))
  g <- default_genotype(m)
  expect_equal(dim(m_matrix(m, g)), c(2L, 2L))
  expect_length(phenotype(m, g), 2)
})

test_that("unstable steady states are flagged", {
  # positive feedback with Hill coefficient 4 is bistable: the middle state is unstable
  m <- grn(c("x1 -> x1"), alpha = 2, theta = 1, hill = 4)
  f <- function(x) 2 * x^4 / (1 + x^4) - x
  mid <- uniroot(f, c(0.3, 1.5))$root
  d <- develop(m, 2, x_start = mid, tmax = 0)  # tmax = 0 still polishes with Newton
  expect_true(d$status %in% c(1L, 2L))
  low <- develop(m, 2, x_start = 0.1)
  high <- develop(m, 2, x_start = 3)
  expect_equal(low$status, 0L)
  expect_equal(high$status, 0L)
  expect_lt(low$X[1, 1], mid)
  expect_gt(high$X[1, 1], mid)
})

test_that("linear maps develop, invert and give M = Sigma B B'", {
  B <- matrix(c(5, -4, -4, 5) / 3, 2)
  m <- gp_linear(B)
  expect_equal(phenotype(m, c(1, 2)), c(x1 = -1, x2 = 2), tolerance = 1e-12)
  expect_equal(genotype_at(m, c(-1, 2)), c(1, 2), tolerance = 1e-12)
  expect_equal(unname(m_matrix(m, c(1, 2), mu = 0.01, sigma = 0.05)),
               0.01 * 0.05^2 * B %*% t(B), tolerance = 1e-12)
})

test_that("Lande separatrix of two symmetric peaks is the midline when G is isotropic", {
  land <- landscape_gaussian(rbind(c(3, 5), c(7, 5)), omega = 2)
  sep <- lande_separatrices(paper_motif(6), land, lower = 0, upper = 10, step = 0.02)
  expect_equal(nrow(attr(sep, "saddles")), 1)
  expect_equal(attr(sep, "saddles")[1, ], c(5, 5), tolerance = 1e-3)
  expect_lt(max(abs(sep$x1 - 5)), 0.02)
  # positive covariance rotates the separatrix towards the leading eigenvector (1, 1)
  B <- matrix(c(5, 4, 4, 5) / 3, 2)
  sep2 <- lande_separatrices(gp_linear(B), land, lower = 0, upper = 10, step = 0.02)
  tilt <- coef(lm(x1 ~ x2, data = sep2))[2]
  expect_gt(tilt, 0.3)
})

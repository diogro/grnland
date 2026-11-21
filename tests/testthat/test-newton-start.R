test_that("the Newton start agrees with integration", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  set.seed(1)
  g0 <- c(genotype_at(paper_motif(3), c(6, 2)), 5)
  G <- t(g0 + matrix(rnorm(3 * 200, 0, 0.2), 3))
  X0 <- develop(m, matrix(g0, 1))$X[rep(1, 200), ]       # parental state
  a <- develop(m, G, X0, newton_radius = 0)
  b <- develop(m, G, X0)
  expect_identical(a$status, b$status)
  expect_equal(b$X, a$X, tolerance = 1e-8)
})

test_that("the Newton start stays on the attractor that integration reaches", {
  # positive feedback with Hill coefficient 4 is bistable for these alphas
  m <- grn(c("x1 -> x1"), alpha = 2, theta = 1, hill = 4, evolvable = "alpha")
  grid <- expand.grid(alpha = seq(1.2, 3, by = 0.02), x0 = seq(0.05, 3, by = 0.05))
  a <- develop(m, matrix(grid$alpha), matrix(grid$x0), newton_radius = 0)
  for (r in c(0.1, 5)) {
    b <- develop(m, matrix(grid$alpha), matrix(grid$x0), newton_radius = r)
    expect_identical(a$status, b$status)
    expect_equal(b$X, a$X, tolerance = 1e-8)
  }
})

test_that("the Newton start never accepts an unstable steady state", {
  m <- grn(c("x1 -> x1"), alpha = 2, theta = 1, hill = 4)
  mid <- uniroot(function(x) 2 * x^4 / (1 + x^4) - x, c(0.3, 1.5))$root
  d <- develop(m, 2, x_start = mid, newton_radius = 1)
  expect_true(d$status %in% c(0L, 2L))
  if (d$status == 0L) expect_gt(abs(d$X[1, 1] - mid), 0.1)
})

test_that("simulations with and without the Newton start agree", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  land <- landscape_gaussian(c(6, 2), omega = 0.5)
  g0 <- c(genotype_at(paper_motif(3), c(6, 2)), 5)
  run <- function(r) simulate_population(m, land, start_genotype = g0, N = 100, generations = 60,
                                         mu = 0.05, sigma = 0.1, reproduction = "sexual", record_every = 20,
                                         record_genotype = TRUE, seed = 4, dev_args = list(newton_radius = r))
  expect_equal(run(0.1)$trajectory, run(0)$trajectory, tolerance = 1e-8)
})

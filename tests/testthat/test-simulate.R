test_that("a clonal population climbs a single Gaussian peak", {
  m <- paper_motif(6)
  land <- landscape_gaussian(c(4, 4), omega = 4)
  sim <- simulate_population(m, land, start = c(3, 3), N = 200, generations = 3000,
                             sigma = 0.1, stop = stop_at_peaks(matrix(c(4, 4), 1), radius = 0.2),
                             seed = 1)
  expect_equal(sim$outcome, 1L)
  expect_lt(sim$generations, 3000)
  last <- sim$trajectory[nrow(sim$trajectory), ]
  expect_equal(c(last$mean_x1, last$mean_x2), c(4, 4), tolerance = 0.1)
})

test_that("simulations are reproducible with a seed", {
  m <- paper_motif(1)
  land <- landscape_gaussian(c(5, 5), omega = 4)
  a <- simulate_population(m, land, start = c(3, 3), N = 50, generations = 50, seed = 3)
  b <- simulate_population(m, land, start = c(3, 3), N = 50, generations = 50, seed = 3)
  expect_identical(a$trajectory, b$trajectory)
})

test_that("diploid sexual populations run and keep genotypes consistent", {
  m <- paper_motif(3)
  land <- landscape_gaussian(c(5, 5), omega = 4)
  sim <- simulate_population(m, land, start = c(3, 3), N = 100, generations = 100,
                             ploidy = 2, reproduction = "sexual", sigma = 0.1,
                             record_m = TRUE, seed = 2)
  expect_equal(dim(sim$final$alleles), c(100L, 2L, 2L))
  expect_equal(sim$final$genotypes, sim$final$alleles[, , 1] + sim$final$alleles[, , 2],
               ignore_attr = TRUE)
  X <- phenotype(m, sim$final$genotypes, x_start = sim$final$expression)
  expect_equal(unname(X), unname(sim$final$expression), tolerance = 1e-6)
  expect_true(all(c("M_x1_x1", "M_x1_x2", "M_x2_x2") %in% names(sim$trajectory)))
})

test_that("stored phenotypes always match genotypes (no stale phenotypes)", {
  m <- paper_motif(5)
  land <- landscape_gaussian(c(6, 6), omega = 4)
  sim <- simulate_population(m, land, start = c(4, 4), N = 100, generations = 200,
                             sigma = 0.2, seed = 5)
  X <- phenotype(m, sim$final$genotypes, x_start = sim$final$expression)
  expect_equal(unname(X), unname(sim$final$expression), tolerance = 1e-6)
})

test_that("Perlin landscapes and peak finding work", {
  skip_if_not_installed("ambient")
  land <- landscape_perlin(frequency = 0.1, seed = 1234)
  pk <- find_peaks(land, 0, 15, res = 200)
  expect_gt(nrow(pk), 0)
  expect_equal(land(as.matrix(pk[, 1:2])), pk$fitness)
})

test_that("grid-cached landscapes match the original closely", {
  skip_if_not_installed("ambient")
  land <- landscape_perlin(frequency = 0.5, seed = 1245)
  cached <- landscape_grid(land, 0, 15, step = 0.005)
  X <- matrix(runif(2000, 0.5, 14.5), ncol = 2)
  expect_lt(max(abs(cached(X) - land(X))), 1e-4)
})

test_that("find_peaks merges duplicate maxima on flat tops", {
  flat <- landscape_function(function(X) pmin(exp(-rowSums((X - 5.0125)^2) / 4), 0.9995), dim = 2)
  expect_equal(nrow(find_peaks(flat, 0, 10, res = 400)), 1)
  expect_gt(nrow(find_peaks(flat, 0, 10, res = 400, merge = 0)), 1)
})

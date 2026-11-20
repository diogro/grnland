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

test_that("record_genotype records genotype means and covariances", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  land <- landscape_gaussian(c(5, 3), omega = 1)
  sim <- simulate_population(m, land, start_genotype = c(8, 9.5, 5), N = 200,
                             generations = 50, mu = 0.05, sigma = 0.1, record_every = 25,
                             selection = selection_proportional(), record_genotype = TRUE, seed = 3)
  tr <- sim$trajectory
  expect_true(all(c("gmean_alpha[x1]", "gmean_theta[x2,x1]",
                    "G_alpha[x1]_alpha[x2]", "G_theta[x2,x1]_theta[x2,x1]") %in% names(tr)))
  expect_equal(unlist(tr[1, paste0("gmean_", loci_names(m))]), c(8, 9.5, 5), ignore_attr = TRUE)
  last <- tr[nrow(tr), ]
  G <- sim$final$genotypes[!sim$final$invalid, , drop = FALSE]
  expect_equal(unlist(last[paste0("gmean_", loci_names(m))]), colMeans(G), ignore_attr = TRUE)
  expect_equal(last[["G_alpha[x1]_theta[x2,x1]"]], cov(G)[1, 3])
})

test_that("record_third_moments records third central moments of genotypes", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  land <- landscape_gaussian(c(5, 3), omega = 1)
  sim <- simulate_population(m, land, start_genotype = c(8, 9.5, 5), N = 200,
                             generations = 50, mu = 0.05, sigma = 0.1, record_every = 25,
                             selection = selection_proportional(), record_third_moments = TRUE,
                             seed = 3)
  tr <- sim$trajectory
  ln <- loci_names(m)
  g3 <- grep("^G3_", names(tr), value = TRUE)
  expect_length(g3, 10)
  expect_false(any(grepl("^gmean_", names(tr))))
  expect_equal(unlist(tr[1, g3]), rep(0, 10), ignore_attr = TRUE)
  G <- sim$final$genotypes[!sim$final$invalid, , drop = FALSE]
  D <- sweep(G, 2, colMeans(G))
  last <- tr[nrow(tr), ]
  expect_equal(last[[sprintf("G3_%s_%s_%s", ln[1], ln[1], ln[3])]], mean(D[, 1]^2 * D[, 3]))
  expect_equal(last[[sprintf("G3_%s_%s_%s", ln[3], ln[3], ln[3])]], mean(D[, 3]^3))
  expect_equal(last[[sprintf("G3_%s_%s_%s", ln[1], ln[2], ln[3])]], mean(D[, 1] * D[, 2] * D[, 3]))
})

test_that("develop_unique gives the same result as developing every individual", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  G <- rbind(c(8, 9.5, 5), c(8, 9.5, 5), c(7, 9, 4), c(8, 9.5, 5), c(7, 9, 4))
  X0 <- matrix(1, nrow(G), m$n_genes)
  dev <- function(G, X0) develop(m, G, X0)
  calls <- 0
  counted <- function(G, X0) { calls <<- calls + nrow(G); dev(G, X0) }
  d <- develop_unique(counted)(G, X0)
  expect_identical(d, dev(G, X0))
  expect_equal(calls, 2)
  # different starting states are developed separately with the exact key
  X0[2, 1] <- 2
  calls <- 0
  expect_identical(develop_unique(counted)(G, X0), dev(G, X0))
  expect_equal(calls, 3)
  # the genotype key ignores starting states
  calls <- 0
  develop_unique(counted, dev_key_genotype)(G, X0)
  expect_equal(calls, 2)
})

test_that("dev keys compare values exactly", {
  x <- 0.1 + 0.2
  G <- rbind(c(x, 1), c(0.3, 1))
  X0 <- matrix(0, 2, 2)
  expect_false(dev_key_exact(G, X0)[1] == dev_key_exact(G, X0)[2])
  expect_equal(dev_key_genotype(rbind(c(1, 2), c(1, 2)), X0)[1], dev_key_genotype(rbind(c(1, 2), c(1, 2)), X0)[2])
})

test_that("develop_unique refuses results it cannot expand", {
  dev <- function(G, X0) list(X = G, status = integer(nrow(G)), extra = 1:5)
  G <- rbind(c(1, 1), c(1, 1), c(2, 2), c(3, 3))
  expect_error(develop_unique(dev)(G, G), "cannot expand")
})

test_that("de-duplicated sexual simulations match the full ones exactly", {
  m <- paper_motif(3, evolvable = c("alpha", "theta[x2,x1]"))
  land <- landscape_gaussian(c(6, 2), omega = 0.5)
  g0 <- c(genotype_at(paper_motif(3), c(6, 2)), 5)
  for (ploidy in 1:2) {
    run <- function(dedupe) simulate_population(m, land, start_genotype = g0, N = 100, generations = 60,
                                                mu = 0.05, sigma = 0.1, ploidy = ploidy, reproduction = "sexual",
                                                record_every = 20, record_genotype = TRUE, seed = 4, dedupe = dedupe)
    expect_identical(run("exact"), run("none"))
  }
})

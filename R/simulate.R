#' Selection schemes
#'
#' Each returns a function mapping fitness values to parental sampling
#' weights. Individuals whose development failed (`NA` fitness) get weight 0.
#'
#' @name selection
NULL

#' @describeIn selection Weights \eqn{w - \min(w) + \epsilon}. Used by the
#'   original analysis scripts. Selection is relative to the least fit
#'   individual, so it is very strong when the population has little fitness
#'   variance.
#' @param eps Weight given to the least fit individual.
#' @export
selection_relative_min <- function(eps = 0.001) {
  function(w) {
    out <- w - min(w, na.rm = TRUE) + eps
    out[is.na(out)] <- 0
    out
  }
}

#' @describeIn selection Weights proportional to fitness (Wright-Fisher).
#'   Fitness must be non-negative.
#' @export
selection_proportional <- function() {
  function(w) {
    w[is.na(w)] <- 0
    if (any(w < 0)) stop("negative fitness under proportional selection")
    w
  }
}

#' @describeIn selection Weights \eqn{\exp(s w)}.
#' @param s Selection strength.
#' @export
selection_exponential <- function(s = 1) {
  function(w) {
    out <- exp(s * (w - max(w, na.rm = TRUE)))
    out[is.na(out)] <- 0
    out
  }
}

#' Stop a simulation when the population reaches a peak
#'
#' Checks every `every` generations whether at least `fraction` of the
#' population lies within `radius` of one of the `peaks`. This reproduces the
#' convergence rule of the original scripts with `shape = "box"`,
#' `radius = 0.1`, `fraction = 0.5`, `every = 10`.
#'
#' @param peaks Matrix or data frame of peak positions (trait columns first;
#'   extra columns such as `fitness` are ignored).
#' @param radius Half-width of the box, or radius of the ball.
#' @param fraction Fraction of the population required near the peak.
#' @param every Check interval in generations.
#' @param shape `"box"` or `"ball"`.
#' @export
stop_at_peaks <- function(peaks, radius = 0.1, fraction = 0.5, every = 10,
                          shape = c("box", "ball")) {
  shape <- match.arg(shape)
  P <- as.matrix(peaks)
  force(radius); force(fraction)
  fun <- function(X) {
    P <- P[, seq_len(ncol(X)), drop = FALSE]
    for (k in seq_len(nrow(P))) {
      D <- abs(sweep(X, 2, P[k, ]))
      near <- if (shape == "box") apply(D <= radius, 1, all) else sqrt(rowSums(D^2)) <= radius
      if (mean(near, na.rm = TRUE) >= fraction) return(k)
    }
    NULL
  }
  structure(list(fun = fun, every = every, peaks = P), class = "grn_stop")
}

#' Simulate a population evolving on an adaptive landscape
#'
#' Individual-based simulation with discrete, non-overlapping generations:
#' selection (parents sampled with fitness-dependent weights), reproduction,
#' mutation, then development of new genotypes into phenotypes.
#'
#' Each locus has `ploidy` allele copies; the genotypic value of a locus is
#' the sum of its alleles. Mutations hit each allele copy with probability
#' `mu` per generation and add a normal deviate with sd `sigma`; allelic
#' values are kept within `[lower, upper]`.
#'
#' With `reproduction = "clonal"` offspring copy their parent. With
#' `"sexual"`, each offspring has two parents (a female and a male if
#' `sexes = TRUE`), and inherits at each locus one randomly chosen allele copy
#' from each parent (free recombination; for haploids, the allele comes from
#' either parent).
#'
#' The original analysis scripts correspond to `ploidy = 1`,
#' `reproduction = "clonal"`, `selection = selection_relative_min()`,
#' `mu = 0.01`, `sigma = 0.05`, `N = 1000`.
#'
#' @param model A `grn` model.
#' @param landscape Function from a phenotype matrix to fitness values.
#' @param start Starting phenotype (trait values). The population starts
#'   monomorphic at the genotype found by [genotype_at()].
#' @param start_genotype Alternatively, the starting genotypic values.
#' @param N Population size.
#' @param generations Maximum number of generations.
#' @param mu Per-allele mutation rate per locus.
#' @param sigma Standard deviation of mutational effects (scalar or per locus).
#' @param ploidy 1 or 2.
#' @param reproduction `"clonal"` or `"sexual"`.
#' @param sexes For sexual reproduction, split the population into equal
#'   numbers of females and males.
#' @param selection A selection scheme, see [selection].
#' @param lower,upper Bounds on allelic values.
#' @param invalid What happens to individuals whose development fails
#'   (no stable steady state): `"lethal"` gives them zero weight;
#'   `"reject"` (clonal only) discards the mutation that caused the failure.
#' @param record_every Record population summaries every this many
#'   generations (the first and last generations are always recorded).
#' @param record_m Also record the M-matrix at the mean genotype.
#' @param record_genotype Also record the mean genotypic value of each locus
#'   (`gmean_<locus>`) and the genotypic covariance matrix of the loci
#'   (`G_<locus>_<locus>`, upper triangle), over individuals that develop.
#' @param record_third_moments Also record the third central moments of the
#'   genotypic values, E[(g_i - mean_i)(g_j - mean_j)(g_k - mean_k)] for
#'   i <= j <= k (`G3_<locus>_<locus>_<locus>`; divisor n, not n - 1), over
#'   individuals that develop.
#' @param stop Optional stopping rule, see [stop_at_peaks()].
#' @param seed Optional RNG seed.
#' @param dev_args List of extra arguments for [develop()].
#' @return An object of class `grn_sim`: a list with `trajectory` (data frame
#'   of generation, trait means, phenotypic (co)variances `P_*`, mean fitness
#'   and number of failed developments; plus `M_*` if `record_m`, and
#'   `gmean_*`, `G_*` if `record_genotype`; `G3_*` if `record_third_moments`), `outcome`
#'   (index of the peak reached, or `NA`), `generations` run, and `final`
#'   (alleles, genotypes, phenotypes and fitness of the last generation).
#'   With clonal reproduction and no environmental noise the phenotypic
#'   covariance matrix equals the genetic covariance matrix G.
#' @export
simulate_population <- function(model, landscape, start = NULL, start_genotype = NULL,
                                N = 1000, generations = 10000, mu = 0.01, sigma = 0.05,
                                ploidy = 1, reproduction = c("clonal", "sexual"),
                                sexes = TRUE, selection = selection_relative_min(),
                                lower = 0, upper = Inf, invalid = c("lethal", "reject"),
                                record_every = 10, record_m = FALSE, record_genotype = FALSE,
                                record_third_moments = FALSE,
                                stop = NULL,
                                seed = NULL, dev_args = list()) {
  reproduction <- match.arg(reproduction)
  invalid <- match.arg(invalid)
  if (invalid == "reject" && reproduction != "clonal") stop("invalid = 'reject' needs clonal reproduction")
  if (!ploidy %in% 1:2) stop("ploidy must be 1 or 2")
  if (!is.null(seed)) set.seed(seed)
  L <- nrow(model$loci)
  traits <- model$traits
  tn <- model$genes[traits]
  sigma <- rep_len(sigma, L)
  dev <- function(G, X0) do.call(develop, c(list(model, G, X0), dev_args))

  # ---- initial population
  if (is.null(start_genotype)) {
    if (is.null(start)) stop("give 'start' or 'start_genotype'")
    start_genotype <- genotype_at(model, start)
  }
  g0 <- as.numeric(start_genotype)
  x_init <- model$x0
  if (!is.null(start)) x_init[traits] <- start
  d0 <- dev(matrix(g0, 1), matrix(x_init, 1))
  if (d0$status != 0) stop("starting genotype does not develop (status ", d0$status, ")")
  # alleles: N x L x ploidy, each copy carrying an equal share
  A <- array(rep(g0 / ploidy, each = N), c(N, L, ploidy))
  X <- matrix(d0$X[1, ], N, model$n_genes, byrow = TRUE)
  W <- rep(landscape(X[1, traits, drop = FALSE]), N)
  bad <- rep(FALSE, N)
  genotype <- function(A) if (ploidy == 1) A[, , 1, drop = TRUE] else A[, , 1] + A[, , 2]

  # ---- recording
  nt <- length(traits)
  cov_idx <- which(upper.tri(diag(nt), diag = TRUE), arr.ind = TRUE)
  cov_names <- sprintf("P_%s_%s", tn[cov_idx[, 1]], tn[cov_idx[, 2]])
  ln <- model$loci$name
  g_idx <- which(upper.tri(diag(L), diag = TRUE), arr.ind = TRUE)
  g3_idx <- as.matrix(expand.grid(i = seq_len(L), j = seq_len(L), k = seq_len(L)))
  g3_idx <- g3_idx[g3_idx[, 1] <= g3_idx[, 2] & g3_idx[, 2] <= g3_idx[, 3], , drop = FALSE]
  g3_idx <- g3_idx[order(g3_idx[, 1], g3_idx[, 2], g3_idx[, 3]), , drop = FALSE]
  rec_cols <- c("generation", paste0("mean_", tn), cov_names, "mean_fitness", "n_invalid",
                if (record_m) sub("^P_", "M_", cov_names),
                if (record_genotype) c(paste0("gmean_", ln),
                                       sprintf("G_%s_%s", ln[g_idx[, 1]], ln[g_idx[, 2]])),
                if (record_third_moments) sprintf("G3_%s_%s_%s", ln[g3_idx[, 1]], ln[g3_idx[, 2]],
                                                  ln[g3_idx[, 3]]))
  rec <- matrix(NA_real_, generations %/% record_every + 2, length(rec_cols),
                dimnames = list(NULL, rec_cols))
  nrec <- 0
  record <- function(gen) {
    Xt <- X[!bad, traits, drop = FALSE]
    C <- if (nrow(Xt) > 1) stats::cov(Xt) else matrix(0, nt, nt)
    row <- c(gen, colMeans(Xt), C[cov_idx], mean(W, na.rm = TRUE), sum(bad))
    if (record_m) {
      M <- tryCatch(m_matrix(model, colMeans(matrix(genotype(A), N)[!bad, , drop = FALSE]),
                             mu = mu, sigma = sigma, ploidy = ploidy),
                    error = function(e) matrix(NA, nt, nt))
      row <- c(row, M[cov_idx])
    }
    if (record_genotype) {
      Gt <- matrix(genotype(A), N)[!bad, , drop = FALSE]
      GC <- if (nrow(Gt) > 1) stats::cov(Gt) else matrix(0, L, L)
      row <- c(row, colMeans(Gt), GC[g_idx])
    }
    if (record_third_moments) {
      Gt <- matrix(genotype(A), N)[!bad, , drop = FALSE]
      D <- sweep(Gt, 2, colMeans(Gt))
      row <- c(row, if (nrow(D) > 0)
        colMeans(D[, g3_idx[, 1], drop = FALSE] * D[, g3_idx[, 2], drop = FALSE] *
                 D[, g3_idx[, 3], drop = FALSE])
        else rep(NA_real_, nrow(g3_idx)))
    }
    nrec <<- nrec + 1
    rec[nrec, ] <<- row
  }
  record(0)

  outcome <- NA_integer_
  gen <- 0
  for (gen in seq_len(generations)) {
    w <- selection(W)
    if (!any(w > 0)) stop("all individuals have zero weight at generation ", gen)

    # ---- reproduction
    if (reproduction == "clonal") {
      idx <- sample.int(N, N, replace = TRUE, prob = w)
      A <- A[idx, , , drop = FALSE]
      X <- X[idx, , drop = FALSE]
      W <- W[idx]
      bad <- bad[idx]
      changed <- rep(FALSE, N)
    } else {
      if (sexes) {
        fem <- seq_len(N) <= N %/% 2
        if (!any(w[fem] > 0) || !any(w[!fem] > 0)) stop("one sex has zero total weight at generation ", gen)
        mothers <- which(fem)[sample.int(sum(fem), N, TRUE, prob = w[fem])]
        fathers <- which(!fem)[sample.int(sum(!fem), N, TRUE, prob = w[!fem])]
      } else {
        mothers <- sample.int(N, N, TRUE, prob = w)
        fathers <- sample.int(N, N, TRUE, prob = w)
      }
      gamete <- function(par) {
        if (ploidy == 1) return(A[par, , 1, drop = FALSE][, , 1])
        pick <- matrix(sample.int(2L, N * L, TRUE), N, L)
        matrix(A[cbind(rep(par, L), rep(seq_len(L), each = N), as.vector(pick))], N, L)
      }
      if (ploidy == 1) {
        from_mother <- matrix(stats::runif(N * L) < 0.5, N, L)
        Anew <- ifelse(from_mother, gamete(mothers), gamete(fathers))
        A <- array(Anew, c(N, L, 1))
      } else {
        A <- array(c(gamete(mothers), gamete(fathers)), c(N, L, 2))
      }
      X <- X[mothers, , drop = FALSE]
      changed <- rep(TRUE, N)
    }

    # ---- mutation
    n_sites <- N * L * ploidy
    nmut <- stats::rbinom(1, n_sites, mu)
    if (invalid == "reject") A_before <- A
    if (nmut > 0) {
      pos <- sample.int(n_sites, nmut)
      locus <- ((pos - 1) %/% N) %% L + 1
      A[pos] <- pmin(pmax(A[pos] + stats::rnorm(nmut, 0, sigma[locus]), lower), upper)
      changed[(pos - 1) %% N + 1] <- TRUE
    }

    # ---- development and fitness of changed individuals
    if (any(changed)) {
      ch <- which(changed)
      Gc <- matrix(genotype(A), N)[ch, , drop = FALSE]
      d <- dev(Gc, X[ch, , drop = FALSE])
      ok <- d$status == 0
      if (invalid == "reject" && any(!ok)) {
        A[ch[!ok], , ] <- A_before[ch[!ok], , , drop = FALSE]
        ch <- ch[ok]
        d$X <- d$X[ok, , drop = FALSE]
        ok <- ok[ok]
      }
      X[ch, ] <- d$X
      bad[ch] <- !ok
      W[ch] <- NA
      if (any(ok)) W[ch[ok]] <- landscape(d$X[ok, traits, drop = FALSE])
    }

    if (gen %% record_every == 0) record(gen)
    if (!is.null(stop) && gen %% stop$every == 0) {
      hit <- stop$fun(X[!bad, traits, drop = FALSE])
      if (!is.null(hit)) {
        outcome <- hit
        break
      }
    }
  }
  if (rec[nrec, "generation"] != gen) record(gen)

  structure(list(
    trajectory = as.data.frame(rec[seq_len(nrec), , drop = FALSE]),
    outcome = outcome,
    generations = gen,
    final = list(alleles = A, genotypes = matrix(genotype(A), N, dimnames = list(NULL, model$loci$name)),
                 expression = X, fitness = W, invalid = bad),
    settings = list(N = N, mu = mu, sigma = sigma, ploidy = ploidy,
                    reproduction = reproduction, start_genotype = g0)
  ), class = "grn_sim")
}

#' @export
print.grn_sim <- function(x, ...) {
  last <- x$trajectory[nrow(x$trajectory), ]
  cat(sprintf("<grn_sim> %d generations, N = %d, %s, ploidy %d\n", x$generations,
              x$settings$N, x$settings$reproduction, x$settings$ploidy))
  m <- last[grep("^mean_", names(last))]
  m <- m[names(m) != "mean_fitness"]
  cat("final mean phenotype:", paste(sprintf("%s = %.3f", sub("^mean_", "", names(m)), unlist(m)), collapse = ", "), "\n")
  cat("outcome (peak reached):", x$outcome, "\n")
  invisible(x)
}

#' Nearest peak to a phenotype
#'
#' @param x Phenotype vector or matrix (one row per phenotype).
#' @param peaks Matrix or data frame of peak positions.
#' @export
nearest_peak <- function(x, peaks) {
  P <- as.matrix(peaks)
  X <- matrix(x, ncol = if (is.null(dim(x))) length(x) else ncol(x))
  P <- P[, seq_len(ncol(X)), drop = FALSE]
  apply(X, 1, function(r) which.min(colSums((t(P) - r)^2)))
}

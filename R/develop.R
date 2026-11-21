#' Develop genotypes into phenotypes
#'
#' Integrates the regulatory dynamics from a starting expression state until
#' a steady state is reached. The steady state is checked for stability.
#'
#' @param model A `grn` model.
#' @param genotypes Matrix of genotypic values, one row per individual and one
#'   column per locus (see [loci_names()]), or a single vector.
#' @param x_start Starting expression states (matrix with one row per
#'   individual, or one vector used for all). Defaults to `model$x0`. When
#'   the network has several stable states, this determines which is reached.
#' @param tol Relative tolerance on \eqn{|dx/dt|} at the steady state.
#' @param tmax,maxsteps Limits on developmental time and integration steps.
#' @param check_stability Check that the steady state is stable.
#' @param newton_radius Hill networks: first try to solve for the steady
#'   state with Newton's method from `x_start`, accepting the result only if
#'   it is stable, within this relative distance (max-norm) of the start, and
#'   in the direction the dynamics point at the start; otherwise integrate.
#'   Much faster when starts are close to the steady state, as for offspring
#'   developing from their parent's expression. `0` always integrates. Other
#'   models ignore it.
#' @return A list with `X`, the full gene-expression matrix, and `status`
#'   (0 = stable steady state, 1 = no convergence, 2 = unstable steady state,
#'   3 = non-finite values).
#' @export
develop <- function(model, genotypes, x_start = NULL, tol = 1e-10,
                    tmax = 1e4, maxsteps = 1e5, check_stability = TRUE,
                    newton_radius = 0.1) {
  UseMethod("develop")
}

#' @export
develop.hill_grn <- function(model, genotypes, x_start = NULL, tol = 1e-10,
                             tmax = 1e4, maxsteps = 1e5, check_stability = TRUE,
                             newton_radius = 0.1) {
  G <- as_genotype_matrix(model, genotypes)
  X0 <- start_matrix(model, x_start, nrow(G))
  l <- model$loci
  out <- hill_develop_cpp(hill_base(model), l$code, l$i - 1L,
                          ifelse(is.na(l$j), 0L, l$j - 1L), G, X0,
                          tol, tmax, as.integer(maxsteps), check_stability,
                          as.numeric(newton_radius))
  colnames(out$X) <- model$genes
  out
}

#' @export
develop.custom_grn <- function(model, genotypes, x_start = NULL, tol = 1e-10,
                               tmax = 1e4, maxsteps = 1e5, check_stability = TRUE,
                               newton_radius = 0.1) {
  G <- as_genotype_matrix(model, genotypes)
  X0 <- start_matrix(model, x_start, nrow(G))
  N <- nrow(G)
  X <- matrix(NA_real_, N, model$n_genes, dimnames = list(NULL, model$genes))
  status <- integer(N)
  for (r in seq_len(N)) {
    p <- custom_params(model, G[r, ])
    res <- custom_develop_one(model, p, X0[r, ], tol, tmax, maxsteps, check_stability)
    X[r, ] <- res$x
    status[r] <- res$status
  }
  list(X = X, status = status)
}

#' Phenotype (trait values) of genotypes
#'
#' Convenience wrapper around [develop()] returning the trait columns only.
#' Individuals whose development fails get `NA`.
#' @inheritParams develop
#' @param ... Passed to [develop()].
#' @export
phenotype <- function(model, genotypes, x_start = NULL, ...) {
  d <- develop(model, genotypes, x_start, ...)
  X <- d$X[, model$traits, drop = FALSE]
  X[d$status != 0, ] <- NA
  if (is.null(dim(genotypes))) X[1, ] else X
}

#' Derivatives of the developmental system at a steady state
#'
#' @param model A `grn` model.
#' @param genotype Vector of genotypic values.
#' @param x Steady-state expression of all genes. Computed if omitted.
#' @return A list with `x`, the Jacobian `J` (\eqn{\partial f/\partial x}),
#'   `dfdp` (\eqn{\partial f/\partial} loci) and the sensitivity matrix
#'   `S` \eqn{= \partial x^*/\partial} loci \eqn{= -J^{-1} \partial f/\partial p}
#'   (all genes by all loci).
#' @export
sensitivity <- function(model, genotype, x = NULL) {
  genotype <- as.numeric(genotype)
  if (is.null(x)) {
    d <- develop(model, genotype)
    if (d$status != 0) warning("development status ", d$status)
    x <- d$X[1, ]
  }
  if (inherits(model, "linear_gp")) {
    J <- -diag(model$n_genes)
    D <- model$B
  } else if (inherits(model, "hill_grn")) {
    l <- model$loci
    dv <- hill_derivs_cpp(hill_base(model), l$code, l$i - 1L,
                          ifelse(is.na(l$j), 0L, l$j - 1L), genotype, x)
    J <- dv$J
    D <- dv$dfdp
  } else {
    p <- custom_params(model, genotype)
    J <- custom_jacobian(model, x, p)
    D <- vapply(model$loci$param, function(nm) {
      h <- 1e-6 * max(1, abs(p[[nm]]))
      pp <- p; pm <- p
      pp[[nm]] <- pp[[nm]] + h
      pm[[nm]] <- pm[[nm]] - h
      (model$rhs(x, pp) - model$rhs(x, pm)) / (2 * h)
    }, numeric(model$n_genes))
    D <- matrix(D, model$n_genes)
  }
  dimnames(J) <- list(model$genes, model$genes)
  dimnames(D) <- list(model$genes, model$loci$name)
  S <- -solve(J, D)
  list(x = stats::setNames(as.numeric(x), model$genes), J = J, dfdp = D, S = S)
}

#' Mutational variance-covariance matrix (M-matrix)
#'
#' \deqn{M = S_T \Sigma S_T^\top,} where \eqn{S_T} holds the rows of the
#' sensitivity matrix belonging to the traits and \eqn{\Sigma} is the
#' per-generation mutational variance of genotypic values,
#' `ploidy * mu * sigma^2` per locus. This is the local (linearised) M-matrix
#' at the given genotype.
#'
#' @param model A `grn` model.
#' @param genotype Genotypic values. Either this or `phenotype` is needed.
#' @param phenotype Trait values; the genotype is found with [genotype_at()].
#' @param mu Per-allele mutation rate per locus (scalar or per locus).
#' @param sigma Standard deviation of mutational effects (scalar or per locus).
#' @param ploidy 1 (haploid) or 2 (diploid).
#' @export
m_matrix <- function(model, genotype = NULL, phenotype = NULL, mu = 0.01,
                     sigma = 0.05, ploidy = 1) {
  if (is.null(genotype)) {
    if (is.null(phenotype)) stop("give either genotype or phenotype")
    genotype <- genotype_at(model, phenotype)
  }
  S <- sensitivity(model, genotype)$S[model$traits, , drop = FALSE]
  L <- ncol(S)
  v <- rep_len(ploidy * mu * sigma^2, L)
  S %*% (v * t(S))
}

#' Genotype that develops into a given phenotype
#'
#' Analytical for Hill networks where every gene is a trait and the loci are
#' exactly the production rates `alpha`. Otherwise solved by Newton iteration
#' on the sensitivity matrix, which needs as many loci as traits.
#'
#' @param model A `grn` model.
#' @param phenotype Target trait values.
#' @param start Starting genotype for the Newton iteration.
#' @param tol Tolerance on the phenotype.
#' @export
genotype_at <- function(model, phenotype, start = NULL, tol = 1e-10) {
  phenotype <- as.numeric(phenotype)
  nt <- length(model$traits)
  if (length(phenotype) != nt) stop("phenotype must have ", nt, " values")
  if (inherits(model, "linear_gp") && nt == model$n_genes)
    return(as.numeric(solve(model$B, phenotype)))
  if (inherits(model, "hill_grn") && nt == model$n_genes &&
      setequal(model$loci$param, "alpha") && nrow(model$loci) == model$n_genes) {
    x <- numeric(model$n_genes)
    x[model$traits] <- phenotype
    i <- model$loci$i
    prodH <- vapply(i, function(k) hill_input(model, k, x), numeric(1))
    g <- model$gamma[i] * x[i] / prodH
    if (any(!is.finite(g))) stop("no genotype produces this phenotype (zero regulatory input)")
    return(unname(g))
  }
  L <- nrow(model$loci)
  if (L != nt) stop("numerical inversion needs as many loci as traits (", L, " vs ", nt, ")")
  g <- if (is.null(start)) default_genotype(model) else as.numeric(start)
  xs <- model$x0
  for (it in 1:100) {
    d <- develop(model, g, x_start = xs)
    if (d$status != 0) stop("development failed during inversion (status ", d$status, ")")
    xs <- d$X[1, ]
    r <- xs[model$traits] - phenotype
    if (max(abs(r)) < tol * (1 + max(abs(phenotype)))) return(g)
    S <- sensitivity(model, g, xs)$S[model$traits, , drop = FALSE]
    step <- solve(S, r)
    lam <- 1
    repeat {
      gn <- g - lam * step
      if (all(gn >= 0) || lam < 1e-4) break
      lam <- lam / 2
    }
    g <- gn
  }
  stop("genotype_at did not converge")
}

# ---- internals --------------------------------------------------------------

hill_base <- function(model) {
  list(alpha = unname(model$alpha), gamma = unname(model$gamma),
       theta = unname(model$theta), hill = unname(model$hill),
       regulation = unname(model$regulation))
}

hill_input <- function(model, i, x) {
  p <- 1
  for (j in seq_len(model$n_genes)) {
    r <- model$regulation[i, j]
    if (r == 0) next
    th <- model$theta[i, j]; h <- model$hill[i, j]
    a <- if (x[j] <= 0) 0 else x[j]^h / (th^h + x[j]^h)
    p <- p * if (r > 0) a else 1 - a
  }
  p
}

as_genotype_matrix <- function(model, genotypes) {
  L <- nrow(model$loci)
  G <- if (is.null(dim(genotypes))) matrix(as.numeric(genotypes), nrow = 1) else as.matrix(genotypes)
  if (ncol(G) != L) stop("genotypes must have ", L, " columns (", paste(model$loci$name, collapse = ", "), ")")
  storage.mode(G) <- "double"
  G
}

start_matrix <- function(model, x_start, N) {
  if (is.null(x_start)) x_start <- model$x0
  if (is.null(dim(x_start))) {
    if (length(x_start) != model$n_genes) stop("x_start must have ", model$n_genes, " values")
    return(matrix(as.numeric(x_start), N, model$n_genes, byrow = TRUE))
  }
  X0 <- as.matrix(x_start)
  if (nrow(X0) != N || ncol(X0) != model$n_genes) stop("x_start has wrong dimensions")
  storage.mode(X0) <- "double"
  X0
}

custom_params <- function(model, g) {
  p <- model$params
  p[model$loci$param] <- g
  p
}

custom_jacobian <- function(model, x, p) {
  if (!is.null(model$jacobian)) return(as.matrix(model$jacobian(x, p)))
  n <- length(x)
  J <- matrix(0, n, n)
  for (k in seq_len(n)) {
    h <- 1e-6 * max(1, abs(x[k]))
    xp <- x; xm <- x
    xp[k] <- xp[k] + h
    xm[k] <- xm[k] - h
    J[, k] <- (model$rhs(xp, p) - model$rhs(xm, p)) / (2 * h)
  }
  J
}

# Same algorithm as the C++ solver (Dormand-Prince + Newton polish), in R.
custom_develop_one <- function(model, p, x, tol, tmax, maxsteps, check_stability) {
  f <- function(x) model$rhs(x, p)
  fn <- function(fx, x) max(abs(fx)) / (1 + max(abs(x)))
  integrate <- function(x, thresh, dt) {
    t <- 0
    k1 <- f(x)
    for (s in seq_len(maxsteps)) {
      r <- fn(k1, x)
      if (!is.finite(r)) return(list(x = x, status = 3L, dt = dt))
      if (r < thresh) return(list(x = x, status = 0L, dt = dt))
      if (t > tmax) break
      k2 <- f(x + dt * (k1 / 5))
      k3 <- f(x + dt * (3 / 40 * k1 + 9 / 40 * k2))
      k4 <- f(x + dt * (44 / 45 * k1 - 56 / 15 * k2 + 32 / 9 * k3))
      k5 <- f(x + dt * (19372 / 6561 * k1 - 25360 / 2187 * k2 + 64448 / 6561 * k3 - 212 / 729 * k4))
      k6 <- f(x + dt * (9017 / 3168 * k1 - 355 / 33 * k2 + 46732 / 5247 * k3 + 49 / 176 * k4 - 5103 / 18656 * k5))
      xn <- x + dt * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 2187 / 6784 * k5 + 11 / 84 * k6)
      k7 <- f(xn)
      err <- dt * (71 / 57600 * k1 - 71 / 16695 * k3 + 71 / 1920 * k4 - 17253 / 339200 * k5 + 22 / 525 * k6 - k7 / 40)
      en <- max(abs(err) / (1e-9 + 1e-7 * pmax(abs(x), abs(xn))))
      if (!is.finite(en)) return(list(x = x, status = 3L, dt = dt))
      if (en <= 1) {
        t <- t + dt
        x <- xn
        k1 <- k7
      }
      dt <- dt * min(5, max(0.2, if (en > 0) 0.9 * en^-0.2 else 5))
    }
    list(x = x, status = 1L, dt = dt)
  }
  r <- integrate(x, 1e-5, 0.1)
  if (r$status != 0) return(r)
  x <- r$x
  xn <- x
  ok <- FALSE
  for (it in 1:30) {
    dx <- tryCatch(solve(custom_jacobian(model, xn, p), -f(xn)), error = function(e) NULL)
    if (is.null(dx)) break
    xn <- xn + dx
    fx <- f(xn)
    if (any(!is.finite(fx))) break
    if (fn(fx, xn) < tol) { ok <- TRUE; break }
  }
  if (ok && max(abs(xn - x)) < 1e-3 * (1 + max(abs(x)))) {
    x <- xn
  } else {
    r <- integrate(x, tol, r$dt)
    if (r$status != 0) return(r)
    x <- r$x
  }
  status <- 0L
  if (check_stability && max(Re(eigen(custom_jacobian(model, x, p), only.values = TRUE)$values)) >= 0)
    status <- 2L
  list(x = x, status = status)
}

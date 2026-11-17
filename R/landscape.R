#' Adaptive landscapes
#'
#' A landscape is a function taking a matrix of phenotypes (one row per
#' individual) and returning a vector of fitness values. Any such function can
#' be used directly; these constructors cover the common cases.
#'
#' @name landscapes
NULL

#' @describeIn landscapes Perlin noise (via the ambient package), in 2 or 3
#'   dimensions.
#' @param frequency Perlin frequency; higher values give more rugged
#'   landscapes.
#' @param seed Perlin seed.
#' @param offset Added to phenotypes before evaluating the noise (scalar or
#'   one value per dimension), e.g. to shift the landscape off the lattice.
#' @param transform Applied to the raw noise in \eqn{[-1, 1]}; the default
#'   rescales to \eqn{[0, 1]}.
#' @param dim Number of phenotypic dimensions.
#' @param ... Further arguments to `ambient::gen_perlin()`.
#' @export
landscape_perlin <- function(frequency, seed, offset = 0,
                             transform = function(z) (z + 1) / 2, dim = 2, ...) {
  if (!requireNamespace("ambient", quietly = TRUE)) stop("landscape_perlin() needs the ambient package")
  if (!dim %in% 2:3) stop("Perlin landscapes support 2 or 3 dimensions")
  offset <- rep_len(offset, dim)
  extra <- list(...)
  f <- function(X) {
    X <- matrix(X, ncol = dim)
    args <- c(list(x = X[, 1] + offset[1], y = X[, 2] + offset[2]),
              if (dim == 3) list(z = X[, 3] + offset[3]),
              list(frequency = frequency, seed = seed), extra)
    transform(do.call(ambient::gen_perlin, args))
  }
  new_landscape(f, dim, list(type = "perlin", frequency = frequency, seed = seed, offset = offset))
}

#' @describeIn landscapes Sum (or maximum) of Gaussian peaks,
#'   \eqn{W(x) = \sum_k h_k \exp(-\frac12 (x-\theta_k)^\top \Omega^{-1} (x-\theta_k))}.
#' @param optima Matrix of peak positions, one row per peak.
#' @param omega Width of the peaks: scalar variance, vector of variances, or
#'   a covariance matrix (shared by all peaks), or a list of matrices.
#' @param heights Peak heights.
#' @param combine `"sum"` or `"max"` over peaks.
#' @export
landscape_gaussian <- function(optima, omega = 1, heights = 1, combine = c("sum", "max")) {
  combine <- match.arg(combine)
  optima <- matrix(optima, ncol = if (is.null(dim(optima))) length(optima) else ncol(optima))
  K <- nrow(optima)
  d <- ncol(optima)
  to_mat <- function(o) if (length(o) == 1) diag(o, d) else if (is.null(dim(o))) diag(o) else o
  Om <- if (is.list(omega)) lapply(omega, to_mat) else rep(list(to_mat(omega)), K)
  Oi <- lapply(Om, solve)
  heights <- rep_len(heights, K)
  f <- function(X) {
    X <- matrix(X, ncol = d)
    v <- vapply(seq_len(K), function(k) {
      D <- sweep(X, 2, optima[k, ])
      heights[k] * exp(-0.5 * rowSums((D %*% Oi[[k]]) * D))
    }, numeric(nrow(X)))
    v <- matrix(v, nrow(X))
    if (combine == "sum") rowSums(v) else apply(v, 1, max)
  }
  new_landscape(f, d, list(type = "gaussian", optima = optima))
}

#' @describeIn landscapes Smooth maximum of several landscapes,
#'   \eqn{\sum_k w_k e^{\beta w_k} / \sum_k e^{\beta w_k}}.
#' @param beta Sharpness; large values approach the maximum.
#' @export
landscape_smoothmax <- function(..., beta = 100) {
  ls <- list(...)
  d <- attr(ls[[1]], "n_dim")
  f <- function(X) {
    W <- vapply(ls, function(l) l(X), numeric(NROW(X)))
    W <- matrix(W, NROW(X))
    E <- exp(beta * (W - apply(W, 1, max)))
    rowSums(E * W) / rowSums(E)
  }
  new_landscape(f, d, list(type = "smoothmax", beta = beta))
}

#' @describeIn landscapes Wrap any function of a phenotype matrix.
#' @param f Function taking a phenotype matrix and returning fitness values.
#' @export
landscape_function <- function(f, dim) new_landscape(f, dim, list(type = "function"))

#' @describeIn landscapes Precompute a 2D landscape on a grid and interpolate
#'   bilinearly. Much faster for landscapes that are slow to evaluate (e.g.
#'   Perlin noise, where each call has a large fixed cost), at the price of a
#'   small interpolation error. Outside the grid, values are clamped to the
#'   nearest edge.
#' @param landscape A 2D landscape to cache.
#' @param lower,upper Grid bounds (scalar or per dimension).
#' @param step Grid spacing.
#' @export
landscape_grid <- function(landscape, lower = 0, upper = 15, step = 0.005) {
  lower <- rep_len(lower, 2)
  upper <- rep_len(upper, 2)
  xs <- seq(lower[1], upper[1], by = step)
  ys <- seq(lower[2], upper[2], by = step)
  nx <- length(xs); ny <- length(ys)
  Z <- matrix(landscape(as.matrix(expand.grid(xs, ys))), nx, ny)
  f <- function(X) {
    X <- matrix(X, ncol = 2)
    u <- pmin(pmax((X[, 1] - lower[1]) / step, 0), nx - 1)
    v <- pmin(pmax((X[, 2] - lower[2]) / step, 0), ny - 1)
    i <- pmin(floor(u), nx - 2); j <- pmin(floor(v), ny - 2)
    wu <- u - i; wv <- v - j
    i <- i + 1; j <- j + 1
    (1 - wu) * (1 - wv) * Z[cbind(i, j)] + wu * (1 - wv) * Z[cbind(i + 1, j)] +
      (1 - wu) * wv * Z[cbind(i, j + 1)] + wu * wv * Z[cbind(i + 1, j + 1)]
  }
  new_landscape(f, 2, list(type = "grid", source = attr(landscape, "info"), step = step))
}

new_landscape <- function(f, dim, info) {
  structure(f, class = c("grn_landscape", "function"), n_dim = dim, info = info)
}

#' @export
print.grn_landscape <- function(x, ...) {
  info <- attr(x, "info")
  cat(sprintf("<landscape: %s, %d dimensions>\n", info$type, attr(x, "n_dim")))
  invisible(x)
}

#' Find local optima of a landscape
#'
#' In two dimensions, evaluates the landscape on a regular grid and returns
#' cells that are at least as high as all eight neighbours, dropping the grid
#' border (as in the original analysis scripts). In other dimensions, runs
#' `optim()` from random starts and merges nearby solutions.
#'
#' @param landscape A landscape function.
#' @param lower,upper Bounds of the search region (scalar or per dimension).
#' @param res Grid resolution per dimension (2D).
#' @param dim Number of dimensions (taken from the landscape if possible).
#' @param n_starts Number of random starts (other dimensions).
#' @return Data frame with peak coordinates and `fitness`.
#' @export
find_peaks <- function(landscape, lower = 0, upper = 15, res = 500,
                       dim = attr(landscape, "n_dim"), n_starts = 200) {
  if (is.null(dim)) stop("give 'dim'")
  lower <- rep_len(lower, dim)
  upper <- rep_len(upper, dim)
  if (dim == 2) {
    xs <- seq(lower[1], upper[1], length.out = res)
    ys <- seq(lower[2], upper[2], length.out = res)
    Z <- matrix(landscape(as.matrix(expand.grid(xs, ys))), res, res)
    P <- matrix(-Inf, res + 2, res + 2)
    P[2:(res + 1), 2:(res + 1)] <- Z
    is_max <- matrix(TRUE, res, res)
    for (di in -1:1) for (dj in -1:1) {
      if (di == 0 && dj == 0) next
      is_max <- is_max & Z >= P[(2:(res + 1)) + di, (2:(res + 1)) + dj]
    }
    is_max[c(1, res), ] <- FALSE
    is_max[, c(1, res)] <- FALSE
    idx <- which(is_max, arr.ind = TRUE)
    return(data.frame(x1 = xs[idx[, 1]], x2 = ys[idx[, 2]], fitness = Z[idx]))
  }
  starts <- matrix(stats::runif(n_starts * dim, lower, upper), ncol = dim, byrow = TRUE)
  sols <- t(apply(starts, 1, function(s) {
    o <- stats::optim(s, function(x) -landscape(matrix(x, 1)), method = "L-BFGS-B",
                      lower = lower, upper = upper)
    c(o$par, -o$value)
  }))
  inside <- apply(sols[, 1:dim, drop = FALSE], 1, function(x) all(x > lower + 1e-6 & x < upper - 1e-6))
  sols <- sols[inside, , drop = FALSE]
  keep <- list()
  scale <- max(upper - lower) * 1e-3
  for (r in seq_len(nrow(sols))) {
    if (!any(vapply(keep, function(k) sqrt(sum((k[1:dim] - sols[r, 1:dim])^2)) < scale, logical(1))))
      keep[[length(keep) + 1]] <- sols[r, ]
  }
  out <- as.data.frame(do.call(rbind, keep))
  names(out) <- c(paste0("x", seq_len(dim)), "fitness")
  out
}

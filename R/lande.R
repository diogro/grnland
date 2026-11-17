#' Separatrices of the Lande-equation flow
#'
#' Treats evolution of the mean phenotype as the flow
#' \eqn{\dot{\bar z} = G(\bar z) \nabla W(\bar z)}, with the G-matrix taken
#' proportional to the local M-matrix of the genotype-phenotype map, and traces
#' the separatrices between the basins of attraction of different optima: the
#' stable manifolds of the saddle points of the landscape. Because G is
#' positive definite, critical points and their types are those of the
#' landscape; only the separatrices move. Two dimensions only.
#'
#' The landscape is evaluated once on a grid; gradient and Hessian come from
#' finite differences on that grid and are interpolated bilinearly (as in the
#' original analysis scripts). Only the direction of the flow matters, so G
#' can be given up to a scalar factor, and \eqn{\nabla \ln W} and
#' \eqn{\nabla W} give the same separatrices.
#'
#' @param model A `grn` model with two traits (used for the default `G`).
#' @param landscape A 2D landscape function.
#' @param lower,upper Bounds of the region (scalar or per dimension).
#' @param step Grid spacing for the landscape derivatives.
#' @param G Function of a phenotype returning the 2x2 G-matrix. Defaults to
#'   the model's M-matrix at that phenotype.
#' @param seed_step Spacing of the starting points used to find critical points.
#' @param ds,max_steps Arc-length step and maximum number of steps when tracing.
#' @return Data frame of separatrix paths with columns `x1`, `x2`, `saddle`,
#'   `branch` and `step`. The saddle points are in `attr(, "saddles")`.
#' @export
lande_separatrices <- function(model, landscape, lower = 0, upper = 12, step = 0.01,
                               G = NULL, seed_step = 0.5, ds = 0.01, max_steps = 4000) {
  lower <- rep_len(lower, 2)
  upper <- rep_len(upper, 2)
  if (is.null(G)) G <- function(x) m_matrix(model, phenotype = x)
  xs <- seq(lower[1], upper[1], by = step)
  ys <- seq(lower[2], upper[2], by = step)
  nx <- length(xs); ny <- length(ys)
  Z <- matrix(landscape(as.matrix(expand.grid(xs, ys))), nx, ny)

  d_rows <- function(M) rbind(M[2, ] - M[1, ], (M[3:nrow(M), ] - M[1:(nrow(M) - 2), ]) / 2,
                              M[nrow(M), ] - M[nrow(M) - 1, ]) / step
  d_cols <- function(M) t(d_rows(t(M)))
  gx <- d_rows(Z); gy <- d_cols(Z)
  hxx <- d_rows(gx); hyy <- d_cols(gy); hxy <- d_rows(gy)

  interp <- function(M, z) {
    u <- min(max((z[1] - lower[1]) / step, 0), nx - 1)
    v <- min(max((z[2] - lower[2]) / step, 0), ny - 1)
    i <- min(floor(u), nx - 2); j <- min(floor(v), ny - 2)
    wu <- u - i; wv <- v - j
    (1 - wu) * (1 - wv) * M[i + 1, j + 1] + wu * (1 - wv) * M[i + 2, j + 1] +
      (1 - wu) * wv * M[i + 1, j + 2] + wu * wv * M[i + 2, j + 2]
  }
  grad <- function(z) c(interp(gx, z), interp(gy, z))
  hess <- function(z) {
    o <- interp(hxy, z)
    matrix(c(interp(hxx, z), o, o, interp(hyy, z)), 2)
  }
  inside <- function(z) all(z >= lower & z <= upper)

  # critical points by Newton iteration from a grid of seeds
  seeds <- as.matrix(expand.grid(seq(lower[1], upper[1], by = seed_step),
                                 seq(lower[2], upper[2], by = seed_step)))
  crit <- list()
  for (r in seq_len(nrow(seeds))) {
    z <- seeds[r, ]
    for (it in 1:50) {
      g <- grad(z)
      if (sqrt(sum(g^2)) < 1e-6) break
      H <- hess(z)
      if (abs(det(H)) < 1e-6) { z <- NULL; break }
      z <- z - solve(H, g)
      if (!inside(z)) { z <- NULL; break }
    }
    if (is.null(z) || sqrt(sum(grad(z)^2)) > 1e-6) next
    if (!any(vapply(crit, function(c) sqrt(sum((c - z)^2)) < 10 * step, logical(1))))
      crit[[length(crit) + 1]] <- z
  }
  # Keep nondegenerate saddles away from the domain edge (where the one-sided
  # differences are unreliable). Degenerate critical points occur on flat
  # ridges, e.g. along lattice mid-lines of Perlin noise, where the landscape
  # is exactly constant in one direction.
  margin <- 5 * step
  is_saddle <- function(z) {
    if (any(z < lower + margin | z > upper - margin)) return(FALSE)
    ev <- eigen(hess(z), symmetric = TRUE, only.values = TRUE)$values
    ev[1] > 0 && ev[2] < 0 && min(abs(ev)) > 1e-2 * max(abs(ev))
  }
  saddles <- lapply(Filter(is_saddle, crit), unname)

  flow <- function(z) {
    if (!inside(z)) return(NULL)
    Gm <- tryCatch(G(z), error = function(e) NULL)
    if (is.null(Gm) || any(!is.finite(Gm))) return(NULL)
    v <- as.numeric(Gm %*% grad(z))
    nv <- sqrt(sum(v^2))
    if (nv < 1e-14) return(NULL)
    v / nv
  }
  trace <- function(z) {
    path <- matrix(NA_real_, max_steps, 2)
    for (i in seq_len(max_steps)) {
      path[i, ] <- z
      k1 <- flow(z); if (is.null(k1)) break
      k2 <- flow(z - 0.5 * ds * k1); if (is.null(k2)) break
      k3 <- flow(z - 0.5 * ds * k2); if (is.null(k3)) break
      k4 <- flow(z - ds * k3); if (is.null(k4)) break
      z <- z - ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6   # backward in time
      if (!inside(z)) break
    }
    path[stats::complete.cases(path), , drop = FALSE]
  }

  out <- list()
  b <- 0
  for (s in seq_along(saddles)) {
    z <- saddles[[s]]
    Gm <- tryCatch(G(z), error = function(e) NULL)
    if (is.null(Gm)) next
    e <- eigen(Gm %*% hess(z))
    stable <- which(Re(e$values) < 0)
    if (length(stable) != 1) next
    v <- Re(e$vectors[, stable])
    v <- v / sqrt(sum(v^2))
    for (sgn in c(-1, 1)) {
      p <- trace(z + sgn * 1e-6 * v)
      if (nrow(p) == 0) next
      b <- b + 1
      out[[b]] <- data.frame(x1 = p[, 1], x2 = p[, 2], saddle = s, branch = b,
                             step = seq_len(nrow(p)))
    }
  }
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(x1 = numeric(), x2 = numeric(), saddle = integer(), branch = integer(), step = integer())
  attr(res, "saddles") <- if (length(saddles)) do.call(rbind, saddles) else matrix(numeric(), 0, 2)
  res
}

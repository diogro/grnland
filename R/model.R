#' Hill-type gene regulatory network
#'
#' Builds a network of `n` genes with dynamics
#' \deqn{dx_i/dt = \alpha_i \prod_j H_{ij}(x_j) - \gamma_i x_i,}
#' where \eqn{H_{ij}(x) = x^h/(\theta^h + x^h)} if gene j activates gene i,
#' \eqn{1 - x^h/(\theta^h + x^h)} if j represses i, and 1 otherwise (with
#' \eqn{\theta = \theta_{ij}}, \eqn{h = h_{ij}}). The phenotype is the steady
#' state reached by development.
#'
#' @param regulation Either an `n x n` matrix with entries -1, 0 or 1,
#'   where `regulation[i, j]` is the effect of gene `j` on gene `i`
#'   (+1 activation, -1 repression), or a character vector of edges such as
#'   `c("x1 -> x2", "x2 -| x1")`.
#' @param genes Gene names, or the number of genes. Defaults to the row names
#'   of `regulation`, or the genes named in the edges (in order of
#'   appearance).
#' @param alpha,gamma Production and decay rates (scalar or length `n`).
#'   `alpha` values are only defaults when alpha is evolvable.
#' @param theta,hill Regulatory thresholds and Hill coefficients (scalar or
#'   `n x n` matrix, indexed like `regulation`).
#' @param evolvable Which parameters are genetically variable loci. Any of
#'   `"alpha"`, `"gamma"`, `"theta"`, `"hill"` (all entries of that kind; for
#'   theta and hill only existing edges), or single entries such as
#'   `"alpha[2]"`, `"theta[2,1]"`, mixed freely.
#' @param traits Genes whose expression is the phenotype under selection
#'   (names or indices). Defaults to all genes.
#' @param x0 Default initial expression state used when no other starting
#'   point is given.
#' @return An object of class `hill_grn`.
#' @examples
#' # motif 3 of the paper: x2 represses x1, x1 activates x2
#' m <- grn(c("x1 -> x2", "x2 -| x1"), theta = 5)
#' # a three-gene repressilator-like cascade with evolvable thresholds
#' m3 <- grn(c("a -| b", "b -| c", "c -> a"), evolvable = c("alpha", "theta"))
#' @export
grn <- function(regulation, genes = NULL, alpha = 1, gamma = 1, theta = 1,
                hill = 1, evolvable = "alpha", traits = NULL, x0 = 1) {
  reg <- parse_regulation(regulation, genes)
  n <- nrow(reg)
  gene_names <- rownames(reg)
  vec <- function(v, what) {
    if (length(v) == 1) v <- rep(v, n)
    if (length(v) != n) stop(what, " must have length 1 or ", n)
    stats::setNames(as.numeric(v), gene_names)
  }
  mat <- function(v, what) {
    if (length(v) == 1) v <- matrix(v, n, n)
    if (!is.matrix(v) || any(dim(v) != n)) stop(what, " must be a scalar or ", n, "x", n, " matrix")
    storage.mode(v) <- "double"
    dimnames(v) <- list(gene_names, gene_names)
    v
  }
  model <- list(
    n_genes = n,
    genes = gene_names,
    regulation = reg,
    alpha = vec(alpha, "alpha"),
    gamma = vec(gamma, "gamma"),
    theta = mat(theta, "theta"),
    hill = mat(hill, "hill"),
    x0 = vec(x0, "x0")
  )
  model$loci <- hill_loci(evolvable, reg)
  model$traits <- resolve_traits(traits, gene_names)
  class(model) <- c("hill_grn", "grn")
  model
}

parse_regulation <- function(regulation, genes) {
  if (is.character(regulation)) {
    m <- regmatches(regulation, regexec("^\\s*(\\S+)\\s*(->|-\\|)\\s*(\\S+)\\s*$", regulation))
    bad <- lengths(m) == 0
    if (any(bad)) stop("cannot parse edges: ", paste(regulation[bad], collapse = ", "),
                       " (use 'a -> b' for activation, 'a -| b' for repression)")
    from <- vapply(m, `[`, "", 2)
    type <- vapply(m, `[`, "", 3)
    to <- vapply(m, `[`, "", 4)
    if (is.null(genes)) genes <- unique(c(rbind(from, to)))
    if (is.numeric(genes) && length(genes) == 1) stop("with edge notation, give gene names, not a count")
    if (!all(c(from, to) %in% genes)) stop("edges name genes not in 'genes'")
    n <- length(genes)
    reg <- matrix(0L, n, n, dimnames = list(genes, genes))
    for (k in seq_along(from)) reg[to[k], from[k]] <- if (type[k] == "->") 1L else -1L
    return(reg)
  }
  if (is.numeric(regulation) && length(regulation) == 1 && is.null(genes)) {
    # number of genes, no interactions
    regulation <- matrix(0L, regulation, regulation)
  }
  reg <- as.matrix(regulation)
  if (nrow(reg) != ncol(reg)) stop("regulation matrix must be square")
  if (!all(reg %in% c(-1, 0, 1))) stop("regulation entries must be -1, 0 or 1")
  storage.mode(reg) <- "integer"
  n <- nrow(reg)
  if (is.null(genes)) genes <- rownames(reg)
  if (is.null(genes) || (is.numeric(genes) && length(genes) == 1)) genes <- paste0("x", seq_len(n))
  if (length(genes) != n) stop("'genes' has wrong length")
  dimnames(reg) <- list(genes, genes)
  reg
}

resolve_traits <- function(traits, genes) {
  if (is.null(traits)) return(seq_along(genes))
  if (is.character(traits)) {
    idx <- match(traits, genes)
    if (anyNA(idx)) stop("unknown trait genes: ", paste(traits[is.na(idx)], collapse = ", "))
    return(idx)
  }
  as.integer(traits)
}

hill_loci <- function(evolvable, reg) {
  n <- nrow(reg)
  genes <- rownames(reg)
  codes <- c(alpha = 1L, gamma = 2L, theta = 3L, hill = 4L)
  out <- list()
  for (e in evolvable) {
    if (e %in% c("alpha", "gamma")) {
      out[[length(out) + 1]] <- data.frame(param = e, i = seq_len(n), j = NA_integer_)
    } else if (e %in% c("theta", "hill")) {
      ij <- which(reg != 0, arr.ind = TRUE)
      ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
      out[[length(out) + 1]] <- data.frame(param = rep(e, nrow(ij)), i = ij[, 1], j = ij[, 2])
    } else {
      m <- regmatches(e, regexec("^(alpha|gamma|theta|hill)\\[\\s*([^],]+)\\s*(?:,\\s*([^]]+)\\s*)?\\]$", e, perl = TRUE))[[1]]
      if (length(m) == 0) stop("cannot parse evolvable parameter '", e, "'")
      idx <- function(s) if (grepl("^[0-9]+$", s)) as.integer(s) else match(s, genes)
      i <- idx(m[3])
      j <- if (m[4] == "") NA_integer_ else idx(m[4])
      if (is.na(i) || i > n) stop("bad index in '", e, "'")
      if (m[2] %in% c("theta", "hill")) {
        if (is.na(j) || j > n) stop("'", e, "' needs two indices")
        if (reg[i, j] == 0) warning("'", e, "' refers to a non-existent edge; it has no effect")
      }
      out[[length(out) + 1]] <- data.frame(param = m[2], i = i, j = j)
    }
  }
  loci <- unique(do.call(rbind, out))
  rownames(loci) <- NULL
  loci$code <- codes[loci$param]
  loci$name <- ifelse(is.na(loci$j),
                      sprintf("%s[%s]", loci$param, genes[loci$i]),
                      sprintf("%s[%s,%s]", loci$param, genes[loci$i], genes[pmax(loci$j, 1)]))
  loci
}

#' The six two-gene motifs of Netz, Melo & Milocco
#'
#' Motifs from Gjuvsland et al. (2007) and Hether & Hohenlohe (2014), with
#' evolvable production rates `alpha`:
#' 1: x1 -> x2; 2: x1 -| x2; 3: x1 -> x2, x2 -| x1; 4: x1 -| x2, x2 -| x1;
#' 5: x1 -> x2, x2 -> x1; 6: no interactions.
#'
#' @param k Motif number, 1 to 6.
#' @param theta,gamma Threshold and decay rate (paper: 5 and 1).
#' @param ... Passed to [grn()], e.g. `evolvable`.
#' @export
paper_motif <- function(k, theta = 5, gamma = 1, ...) {
  edges <- switch(as.character(k),
    "1" = "x1 -> x2",
    "2" = "x1 -| x2",
    "3" = c("x1 -> x2", "x2 -| x1"),
    "4" = c("x1 -| x2", "x2 -| x1"),
    "5" = c("x1 -> x2", "x2 -> x1"),
    "6" = character(0),
    stop("motif must be 1 to 6"))
  reg <- matrix(0L, 2, 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
  if (length(edges)) reg <- parse_regulation(edges, genes = c("x1", "x2"))
  grn(reg, theta = theta, gamma = gamma, ...)
}

#' User-defined regulatory network
#'
#' For regulatory relations that the Hill model in [grn()] cannot express.
#' Development and derivatives run in R, so this is slower than [grn()].
#'
#' @param rhs Function `rhs(x, p)` returning \eqn{dx/dt} for expression vector
#'   `x` and named parameter vector `p`.
#' @param params Named numeric vector of all parameters (defaults).
#' @param genes Gene names, or number of genes.
#' @param evolvable Names of parameters in `params` that are loci.
#' @param traits Genes under selection (names or indices).
#' @param jacobian Optional function `jacobian(x, p)` returning the
#'   \eqn{n \times n} matrix \eqn{\partial f/\partial x}. Numerical if omitted.
#' @param x0 Default initial expression state.
#' @return An object of class `custom_grn`.
#' @examples
#' # additive (OR-logic) activation of x2 by x1, which grn() cannot express
#' m <- grn_custom(
#'   function(x, p) c(p["a1"] - x[1],
#'                    p["a2"] * (0.2 + x[1] / (5 + x[1])) - x[2]),
#'   params = c(a1 = 2, a2 = 2), genes = c("x1", "x2"))
#' @export
grn_custom <- function(rhs, params, genes, evolvable = names(params),
                       traits = NULL, jacobian = NULL, x0 = 1) {
  if (is.null(names(params))) stop("params must be named")
  if (!all(evolvable %in% names(params))) stop("evolvable must name entries of params")
  if (is.numeric(genes) && length(genes) == 1) genes <- paste0("x", seq_len(genes))
  n <- length(genes)
  if (length(x0) == 1) x0 <- rep(x0, n)
  model <- list(
    n_genes = n,
    genes = genes,
    rhs = function(x, p) as.numeric(rhs(x, p)),
    jacobian = jacobian,
    params = params,
    loci = data.frame(param = evolvable, name = evolvable),
    traits = resolve_traits(traits, genes),
    x0 = stats::setNames(as.numeric(x0), genes)
  )
  class(model) <- c("custom_grn", "grn")
  model
}

#' @export
print.grn <- function(x, ...) {
  cat(sprintf("<%s> %d genes, %d loci, traits: %s\n", class(x)[1], x$n_genes,
              nrow(x$loci), paste(x$genes[x$traits], collapse = ", ")))
  if (inherits(x, "hill_grn")) {
    ij <- which(x$regulation != 0, arr.ind = TRUE)
    if (nrow(ij)) {
      e <- sprintf("%s %s %s", x$genes[ij[, 2]],
                   ifelse(x$regulation[ij] > 0, "->", "-|"), x$genes[ij[, 1]])
      cat("edges:", paste(e, collapse = ", "), "\n")
    } else cat("edges: none\n")
  }
  cat("loci:", paste(x$loci$name, collapse = ", "), "\n")
  invisible(x)
}

#' Names of a model's loci
#' @param model A `grn` model.
#' @export
loci_names <- function(model) model$loci$name

#' Default genotype of a model (current values of the evolvable parameters)
#' @param model A `grn` model.
#' @export
default_genotype <- function(model) {
  if (inherits(model, "custom_grn")) return(unname(model$params[model$loci$param]))
  if (inherits(model, "linear_gp")) return(rep(1, nrow(model$loci)))
  l <- model$loci
  vapply(seq_len(nrow(l)), function(k) {
    switch(l$param[k],
           alpha = model$alpha[l$i[k]], gamma = model$gamma[l$i[k]],
           theta = model$theta[l$i[k], l$j[k]], hill = model$hill[l$i[k], l$j[k]])
  }, numeric(1))
}

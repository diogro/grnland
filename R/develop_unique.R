#' Develop each distinct individual once
#'
#' Wraps a development function so that individuals with the same key are
#' developed only once and share the result. In sexual populations most
#' offspring repeat a genotype that is already present, so this saves most of
#' the development time. In clonal populations only mutants are developed,
#' and they are almost always distinct, so it changes little.
#'
#' The wrapper is independent of the simulation loop. It applies to any
#' function `dev(G, X0)` that returns a list whose elements have one row (or
#' one entry) per individual, as [develop()] does (`X`, `status`). Other
#' kinds of individuals, for example several motifs in one population, need
#' a `key` that tells them apart (e.g. a motif id), or no de-duplication.
#'
#' @param dev A function `dev(G, X0)` of a genotype matrix and a matrix of
#'   starting expression states, one row per individual.
#' @param key A function `key(G, X0)` returning one string per individual.
#'   Individuals with equal keys must develop identically. See [dev_key_exact()]
#'   and [dev_key_genotype()].
#' @return A function with the same interface as `dev`.
#' @seealso [simulate_population()] (argument `dedupe`).
#' @export
develop_unique <- function(dev, key = dev_key_exact) {
  force(dev)
  force(key)
  function(G, X0) {
    k <- key(G, X0)
    first <- !duplicated(k)
    if (all(first)) return(dev(G, X0))
    d <- dev(G[first, , drop = FALSE], X0[first, , drop = FALSE])
    idx <- match(k, k[first])
    nu <- sum(first)
    lapply(d, function(v) {
      if (is.matrix(v) && nrow(v) == nu) v[idx, , drop = FALSE]
      else if (is.null(dim(v)) && length(v) == nu) v[idx]
      else stop("develop_unique: cannot expand an element of the development result")
    })
  }
}

#' Keys for de-duplicated development
#'
#' `dev_key_exact()` keys on the exact genotype and the exact starting
#' expression state, so de-duplicated development gives bit-identical results
#' to developing every individual. `dev_key_genotype()` keys on the genotype
#' only. Each genotype is then developed once from the starting state of its
#' first occurrence. That is also exact when each genotype has a single stable
#' state reachable from all the starting states present (the usual case away
#' from bistable networks), and otherwise can pick a different steady state.
#' Values are compared exactly (hexadecimal representation of the doubles).
#'
#' @param G Genotype matrix, one row per individual.
#' @param X0 Starting expression states, one row per individual.
#' @return A character vector with one key per row.
#' @export
dev_key_exact <- function(G, X0) row_keys(cbind(G, X0))

#' @rdname dev_key_exact
#' @export
dev_key_genotype <- function(G, X0) row_keys(G)

row_keys <- function(M) {
  h <- matrix(sprintf("%a", M), nrow(M))
  do.call(paste, c(lapply(seq_len(ncol(h)), function(j) h[, j]), sep = "|"))
}

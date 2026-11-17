#' Linear genotype-phenotype map
#'
#' Phenotype \eqn{x = B \alpha}, the linear map used in Fig 1 of the paper
#' (Wagner 1984). Its M-matrix is constant: \eqn{M = \Sigma B B^\top} with
#' \eqn{\Sigma} = `ploidy * mu * sigma^2`.
#'
#' @param B Square matrix mapping genotypic values to traits.
#' @param genes Trait names.
#' @return An object of class `linear_gp`, usable wherever a `grn` model is.
#' @examples
#' B2 <- matrix(c(5, 4, 4, 5) / 3, 2)
#' m <- gp_linear(B2)
#' m_matrix(m, c(1, 1), mu = 0.01, sigma = 0.05)
#' @export
gp_linear <- function(B, genes = paste0("x", seq_len(nrow(B)))) {
  B <- as.matrix(B)
  if (nrow(B) != ncol(B)) stop("B must be square")
  n <- nrow(B)
  dimnames(B) <- list(genes, paste0("a", seq_len(n)))
  model <- list(n_genes = n, genes = genes, B = B,
                loci = data.frame(param = colnames(B), name = colnames(B)),
                traits = seq_len(n), x0 = stats::setNames(rep(0, n), genes))
  class(model) <- c("linear_gp", "grn")
  model
}

#' @export
develop.linear_gp <- function(model, genotypes, x_start = NULL, ...) {
  G <- as_genotype_matrix(model, genotypes)
  X <- G %*% t(model$B)
  colnames(X) <- model$genes
  list(X = X, status = integer(nrow(G)))
}

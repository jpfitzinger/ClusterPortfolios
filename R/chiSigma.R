#' @name chiSigma
#' @title A convex hierarchical filter of the covariance matrix.
#' @description Calculates covariance matrix filtered using the most diversified hierarchical graph.
#' The hierarchical graph can be optimized using maximum diversification (MaxDiv) or equal risk contribution (ERC).
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#'
#' The argument \code{meta_loss} represents the loss function used to optimize the most diversified hierarchical allocation graph.
#' The optimized hierarchy is used to filter \code{sigma}. If the filtered covariance matrix is used in a
#' minimum variance portfolio optimizer, a CHI portfolio is constructed.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @param meta_loss a loss function of the most diversified hierarchical allocation graph.
#' @return A \eqn{(N \times N)}{(N x N)} filtered covariance matrix.
#' @author Johann Pfitzinger
#' @references
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' chiSigma(sigma)
#'
#' @export


chiSigma <- function(
  sigma,
  cluster_method = c("single", "average", "complete", "ward", "DIANA"),
  meta_loss = c("MaxDiv", "ERC")
) {

  cluster_method <- match.arg(cluster_method)
  meta_loss <- match.arg(meta_loss)

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  cluster_object <- .get_clusters(sigma, cluster_method)

  lvl_w <- function(n_clusters, sigma, cluster_object) {

    cut <- cutree(cluster_object, k = n_clusters)
    max_cut <- max(cut)

    cut_fx <- function(rowSel,cut) as.data.frame(matrix(as.numeric(rowSel == cut), ncol = length(cut)))
    S_Filler <- purrr::map(1:max_cut, ~cut_fx(., cut))
    S = matrix(nrow = length(S_Filler), ncol = length(cut))
    # Can be improved with Rcpp if required.
    for(i in 1:length(S_Filler) ) S[i,] <- as.matrix(S_Filler[i][[1]])

    S_av <- sweep(S, 1, rowSums(S), "/")
    sigma_sub <- S_av %*% sigma %*% t(S_av)

    Amat <- cbind(1, -diag(max_cut), diag(max_cut))
    bvec <- c(1, -rep(1, max_cut), rep(0, max_cut))

    opt <- quadprog::solve.QP(sigma_sub, rep(0, max_cut), Amat, bvec, meq = 1)
    w <- as.numeric(t(S_av) %*% opt$solution)

    sigma_ <- cbind(rbind(sigma, 1), c(rep(1, nrow(sigma)), 0))
    S_av_ <- as.matrix(Matrix::bdiag(S_av, 1))
    sigma_sub_ <- cbind(rbind(sigma_sub, 1), c(rep(1, nrow(sigma_sub)), 0))
    rot_mat <- t(S_av_) %*% solve(sigma_sub_) %*% S_av_ %*% sigma_

    return(list(w = w, rotation = rot_mat))

  }

  level_k <- c(n:1)
  w_opt_lvl <- purrr::map(level_k, ~lvl_w(., sigma, cluster_object$cluster_object))
  w_mat <- sapply(w_opt_lvl, function(x) x$w)

  sigma_meta <- t(w_mat) %*% sigma %*% w_mat
  diag(sigma_meta) <- diag(sigma_meta) + 1e-5

  if (meta_loss == "MaxDiv") {

    w_bounds <- diag(n)
    w_bounds[upper.tri(w_bounds)] <- 1

    Amat <- cbind(sqrt(diag(sigma_meta)) / mean(sqrt(diag(sigma_meta))), 1 , -w_bounds, w_bounds)
    bvec <- c(1, 0.9999, rep(-0.99999, n), rep(0.00001, n))
    dvec <- rep(0, n)

    meta_opt <- quadprog::solve.QP(sigma_meta,
                                   dvec,
                                   Amat, bvec, meq = 2)

    phi <- meta_opt$solution / sum(meta_opt$solution)

  }

  if (meta_loss == "ERC") {

    .pRC <- function(w, w_mat, sigma) {
      w <- as.numeric(w_mat %*% w)
      sigmaw <- crossprod(sigma, w)
      pRC <- (w * sigmaw)/as.numeric(crossprod(w, sigmaw))
      d <- sum((pRC - 1/n)^2)
      return(d)
    }
    .eqConstraint <- function (w)
    {
      return(sum(w) - 1)
    }
    phi <- nloptr::slsqp(x0 = rep(1/n, n), fn = .pRC,
                         heq = .eqConstraint, lower = rep(1e-5, n),
                         upper = rep(1, n), nl.info = FALSE,
                         control = list(xtol_rel = 1e-18, check_derivatives = FALSE,
                                        maxeval = 20000),
                         sigma = sigma, w_mat = w_mat)$par

  }

  rot_mat <- lapply(w_opt_lvl, function(x) x$rotation)
  rot_mat <- abind::abind(rot_mat, along = 3)
  rot_mat <- apply(rot_mat, 2, function(x) x %*% phi)
  rot_mat <- solve(rot_mat[1:n, 1:n])

  chiSigma_Mat <- t(rot_mat) %*% sigma %*% rot_mat
  diag(chiSigma_Mat) <- diag(chiSigma_Mat) + 1e-8

  return(chiSigma_Mat)

}


#' @name NCO
#' @title Nested Clusters Optimization
#' @description Computes optimal NCO portfolio with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#' The number of clusters can be passed using the \code{n_clusters} argument,
#' calculated automatically with \code{n_clusters='auto'} using the Silhouette criterion.
#'
#' NCO calculates within cluster and between cluster minimum variance portfolios.
#' Constraints are implemented using an iterative convergence algorithm.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param n_clusters trade-off between naive (0) or cluster-based (1) tree-splitting (see Details).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author Johann Pfitzinger
#' @references
#' Lopez de Prado, M. (2019).
#' A Robust Estimator of the Efficient Frontier.
#' \emph{SSRN Electronic Journal}.
#'
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' NCO(sigma, UB = 0.15, n_clusters = 'auto')
#'
#' @export

NCO <- function(
  sigma,
  cluster_method = c("single", "average", "complete", "ward", "DIANA"),
  UB = NULL,
  LB = NULL,
  n_clusters = "auto"
  ) {

  cluster_method <- match.arg(cluster_method)

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  if (!is.integer(n_clusters) & !is.numeric(n_clusters) & n_clusters != "auto")
    stop("'n_clusters' must be an integer!")

  if (n_clusters != "auto" & n_clusters > n)
    stop("'n_clusters' must be smaller than the number of assets.")

  # Fetch constraints
  if (is.null(UB)) {
    UB <- rep(1, n)
  } else if (length(UB) == 1) {
    # Check constraint
    if (UB * n < 1) stop("Inconsistent constraint (increase UB)")
    UB <- rep(UB, n)
  } else {
    # Check constraint
    if (length(UB) != n) stop("Inconsistent contraint (incorrect elements in UB)")
    UB <- UB
  }
  if (is.null(LB)) {
    LB <- rep(0, n)
  } else if (length(LB) == 1) {
    # Check constraint
    if (LB * n > 1) stop("Inconsistent constraint (decrease LB)")
    LB <- rep(LB, n)
  } else {
    # Check constraint
    if (length(LB) != n) stop("Inconsistent contraint (incorrect elements in LB)")
    LB <- LB
  }
  # Check constraint
  if (!all(pmax(UB, LB) == UB) || !all(pmin(UB, LB) == LB))
    stop("Inconsistent constraint (UB smaller than LB)")

  clust <- .get_clusters(sigma, cluster_method, n_clusters = n_clusters)

  cut <- clust$clusters
  max_cut <- max(cut)

  cut_fx <- function(rowSel,cut) as.data.frame(matrix(as.numeric(rowSel == cut), ncol = length(cut)))
  S_Filler <- purrr::map(1:max_cut, ~cut_fx(., cut))
  S = matrix(nrow = length(S_Filler), ncol = length(cut))
  for(i in 1:length(S_Filler) ) S[i,] <- as.matrix(S_Filler[i][[1]])

  # Inner (intra-level optimization)
  intra_weights <- function(w_inter) {

    S_av <- c()

    for (i in 1:nrow(S)) {

      if (sum(S[i,] == 1) > 1) {

        idx <- S[i,] == 1

        sigma_sub <- sigma[idx, idx]
        n_sub <- sum(idx)

        UB_inner <- UB[idx] / min(sum(UB[idx]), max(sum(w_inter[idx]), 0.001))
        LB_inner <- LB[idx] / max(sum(LB[idx]), max(sum(w_inter[idx]), 0.001))

        w <- MV(sigma_sub, UB = UB_inner, LB = LB_inner)

      } else {

        w <- 1

      }

      S_i <- S[i,]
      S_i[S[i,] == 1] <- w
      S_av <- rbind(S_av, S_i)

    }

    return(S_av)

  }

  inter_weights <- function(S_av) {

    sigma_sub <- S_av %*% sigma %*% t(S_av)

    if (all(dim(sigma_sub) == 1)) {

      opt_weights <- as.numeric(S_av)

    } else {

      # Constraints
      Amat <- cbind(1, -diag(max_cut), diag(max_cut))
      bvec <- c(1, S %*% -UB, S %*% LB)

      # Optimization
      opt <- quadprog::solve.QP(sigma_sub, rep(0, max_cut), Amat, bvec, meq = 1)

      opt_weights <- opt$solution
      opt_weights <- as.numeric(t(S_av) %*% opt_weights)

    }

    return(opt_weights)

  }

  chk <- T
  counter <- 0
  maxit <- 100
  w_init <- rep(1, length(cut))
  while (chk) {

    counter <- counter + 1
    S_av <- intra_weights(w_init)
    opt_weights <- inter_weights(S_av)
    chk <- sum(abs(opt_weights - w_init)) > 0.001 && counter < maxit
    w_init <- opt_weights

  }

  names(opt_weights) <- asset_names

  return(opt_weights)

}


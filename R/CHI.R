#' @name CHI
#' @title Convex Hierarchical Portfolio
#' @description Computes optimal CHI portfolio with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#'
#' The argument \code{meta_loss} represents the loss function used to optimize the most diversified hierarchical allocation graph.
#' The optimized hierarchy is used to filter \code{sigma} and \code{mu}. If the filtered covariance matrix is used in a
#' mean variance portfolio optimizer, a CHI portfolio is constructed.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of estimated returns.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy
#' @param meta_loss a loss function of the most diversified hierarchical allocation graph..
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param gamma risk aversion parameter. Default: \code{gamma = 0}.
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author Johann Pfitzinger
#' @references
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' CHI(sigma, UB = 0.15)
#'
#' @export


CHI <- function(
  sigma,
  mu = NULL,
  cluster_method = c("single", "average", "complete", "ward", "DIANA"),
  meta_loss = c("MaxDiv", "ERC"),
  UB = NULL,
  LB = NULL,
  gamma = 0
  ) {

  cluster_method <- match.arg(cluster_method)
  meta_loss <- match.arg(meta_loss)

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  if (!is.null(mu)) {
    if (length(mu)!=n) {
      stop("Different dimensions implied by 'sigma' and 'mu'")
    }
  }

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

  chi <- chiSigma(sigma, mu, cluster_method, meta_loss, gamma)

  w <- MV(chi$sigma, chi$mu, UB, LB, gamma)

  return(w)

}


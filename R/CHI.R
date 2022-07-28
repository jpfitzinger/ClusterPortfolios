#' @name CHI
#' @title Convex Hierarchical Portfolio
#' @description Computes the optimal CHI-MVO portfolio with full investment, weight and group constraints.
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
#' @param meta_loss a loss function of the most diversified hierarchical allocation graph.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraints.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraints.
#' @param groups vector of group IDs. The names of the vector must be identical to the asset names.
#' @param group.UB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of upper bound group constraints.
#' @param group.LB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of lower bound group constraints.
#' @param gamma risk aversion parameter. Default: \code{gamma = 0} returns the minimum variance portfolio.
#' @param max_tilt maximum percentage reduction in the effective number of assets. Default: \code{max_tilt = 1} (no restriction).
#' @param ... arguments passed to \code{cluster::agnes} method.
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
  meta_loss = c("MaxDiv", "ERC"),
  UB = NULL,
  LB = NULL,
  groups = NULL,
  group.UB = NULL,
  group.LB = NULL,
  gamma = 0,
  max_tilt = 1,
  ...
  ) {

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

  chi <- chiSigma(sigma, mu, meta_loss, UB, LB, gamma, max_tilt, ...)
  # w <- MV(sigma = chi$sigma, mu = chi$mu, UB = UB, LB = LB, gamma = gamma,
  #         groups = groups, group.UB = group.UB, group.LB = group.LB)
  w <- MV(sigma = chi$sigma, mu = mu, UB = UB, LB = LB, gamma = gamma,
          groups = groups, group.UB = group.UB, group.LB = group.LB)
  # w <- chi$w

  return(w)

}


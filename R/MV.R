#' @name MV
#' @title Mean Variance
#' @description Computes a Mean Variance portfolio with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' The MV solution is calculated using \code{quadprog}. When \code{gamma} is left at the
#' default setting, the minimum variance portfolio is computed.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of estimated returns.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param gamma risk aversion parameter. Default: \code{gamma = 0}.
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' MV(sigma, UB = 0.15)
#'
#' @export

MV <- function(
  sigma,
  mu = NULL,
  UB = NULL,
  LB = NULL,
  gamma = 0
) {

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

  if (all(dim(sigma) == 1)) {

    opt_weights <- 1

  } else {

    # Constraints
    Amat <- cbind(1, -diag(n), diag(n))
    bvec <- c(1, -UB, LB)

    if (!is.null(mu)) {
      dvec <- mu
    } else {
      dvec <- rep(0, n)
    }

    # Optimization
    opt <- quadprog::solve.QP(sigma, dvec * gamma, Amat, bvec, meq = 1)

    opt_weights <- opt$solution

  }

  names(opt_weights) <- asset_names

  return(opt_weights)

}


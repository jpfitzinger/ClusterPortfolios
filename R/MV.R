#' @name MV
#' @title Minimum Variance
#' @description Computes a Minimum Variance portfolio with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' The MV solution is calculated using \code{quadprog}.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
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
  UB = NULL,
  LB = NULL
) {

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

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

    # Optimization
    opt <- quadprog::solve.QP(sigma, rep(0, n), Amat, bvec, meq = 1)

    opt_weights <- opt$solution

  }

  names(opt_weights) <- asset_names

  return(opt_weights)

}


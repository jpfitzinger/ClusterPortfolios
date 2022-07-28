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
#' @param groups vector of group IDs. The names of the vector must be identical to the asset names.
#' @param group.UB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of upper bound group constraints.
#' @param group.LB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of lower bound group constraints.
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
  groups = NULL,
  group.UB = NULL,
  group.LB = NULL,
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

  if (!is.null(groups)) {

    n_groups <- length(unique(groups))
    if (length(groups) != n) stop("'groups' has incorrect number of elements")
    if (!all(names(groups) %in% asset_names)) stop("group names must be identical to asset names")
    groups <- groups[asset_names]

    # Fetch constraints
    if (is.null(group.UB)) {
      group.UB <- rep(1, n_groups)
      names(group.UB) <- unique(groups)
    } else if (length(group.UB) == 1) {
      # Check constraint
      if (group.UB * n_groups < 1) stop("Inconsistent constraint (increase group.UB)")
      group.UB <- rep(group.UB, n_groups)
      names(group.UB) <- unique(groups)
    } else {
      # Check constraint
      if (length(group.UB) != n_groups) stop("Inconsistent contraint (incorrect elements in group.UB)")
      group.UB <- group.UB
    }
    if (is.null(group.LB)) {
      group.LB <- rep(0, n_groups)
      names(group.LB) <- unique(groups)
    } else if (length(group.LB) == 1) {
      # Check constraint
      if (group.LB * n_groups > 1) stop("Inconsistent constraint (decrease group.LB)")
      group.LB <- rep(group.LB, n_groups)
      names(group.LB) <- unique(groups)
    } else {
      # Check constraint
      if (length(group.LB) != n_groups) stop("Inconsistent contraint (incorrect elements in group.LB)")
      group.LB <- group.LB
    }

    if (!all(groups %in% names(group.UB)) | !all(groups %in% names(group.LB)))
      stop("Inconsistent constraint (missing group names in 'group.UB' or 'group.LB')")
    group.UB <- group.UB[unique(groups)]
    group.LB <- group.LB[unique(groups)]
    if (!all(pmax(group.UB, group.LB) == group.UB) || !all(pmin(group.UB, group.LB) == group.LB))
      stop("Inconsistent constraint (group.UB smaller than group.LB)")

    groups_mat <- sapply(unique(groups), function(x) x==groups)
    groups_mat <- cbind(-groups_mat, groups_mat)
    group.UB <- -group.UB

  } else {
    groups_mat <- NULL
  }

  if (all(dim(sigma) == 1)) {

    opt_weights <- 1

  } else {

    # Constraints
    Amat <- cbind(1, -diag(n), diag(n), groups_mat)
    bvec <- c(1, -UB, LB, group.UB, group.LB)

    if (!is.null(mu)) {
      dvec <- mu
    } else {
      dvec <- rep(0, n)
    }

    # With return target
    safeOpt <- purrr::safely(quadprog::solve.QP)
    Amat <- cbind(1, -dvec, -diag(n), diag(n), groups_mat)
    bvec <- c(1, -gamma, -UB, LB, group.UB, group.LB)
    opt_UB <- safeOpt(sigma, dvec, Amat, bvec, meq = 1)
    Amat <- cbind(1, dvec, -diag(n), diag(n), groups_mat)
    bvec <- c(1, gamma, -UB, LB, group.UB, group.LB)
    opt_LB <- safeOpt(sigma, -dvec, Amat, bvec, meq = 1)
    if (!is.null(opt_UB$result)) {
      opt <- opt_UB$result
    } else {
      opt <- opt_LB$result
    }


    # Optimization
    # opt <- quadprog::solve.QP(sigma, dvec * gamma, Amat, bvec, meq = 1)

    opt_weights <- opt$solution

  }

  names(opt_weights) <- asset_names

  return(opt_weights)

}


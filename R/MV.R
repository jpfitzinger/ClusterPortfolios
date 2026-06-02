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
#' @author Johann Pfitzinger & Nico Katzke
#' @examples
#' data("spec")
#' spec$n -> n; spec$gamma -> gamma; spec$sigma -> sigma; spec$mu -> mu; spec$group.LB->group.LB; spec$group.UB -> group.UB;
#' spec$groups_mat -> groups_mat; spec$groups -> groups; spec$LB -> LB; spec$UB -> UB; spec$assets -> assets
#' MV( sigma, mu = mu, UB = UB, LB = LB, groups = groups, group.UB = group.UB, group.LB = group.LB, groups_mat = groups_mat, gamma = gamma)
#' MV_Roll_gamma <- seq(0.1, 150, length.out = 100) %>% as.list() %>% map_df(~MV( sigma, mu = mu, UB = UB, LB = LB, groups = groups, group.UB = group.UB, group.LB = group.LB, groups_mat = groups_mat, gamma = .))
#' @export

MV <- function(
    sigma,
    mu = NULL,
    UB = NULL,
    LB = NULL,
    groups = NULL,
    group.UB = NULL,
    group.LB = NULL,
    groups_mat = NULL,
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


  # GROUPS


  if (!is.null(groups)) {

    if(!is.null(groups_mat) & !is.null(colnames(groups_mat))) groups <- colnames(groups_mat)
    n_groups <- length(unique(groups))

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

    if (!all(groups %in% names(group.UB)) | !all(groups %in% names(group.UB))) stop("Inconsistent constraint (missing group names in 'group.UB' or 'group.LB')")

    # Reordering messes with Amat later. remove:
    # group.UB <- group.UB[unique(groups)]
    # group.LB <- group.LB[paste0("LB_", unique(groups))]

    if (!all(pmax(group.UB, group.LB) == group.UB) || !all(pmin(group.UB, group.LB) == group.LB))
      stop("Inconsistent constraint (group.UB smaller than group.LB)")

    #Previous - deprecated, as this does not allow securities to be in different groups at the same time.
    # E.g. a Equity and Global label to an aset class
    #         groups_mat <- sapply(unique(groups), function(x) x==groups)
    #         groups_mat <- cbind(-groups_mat, groups_mat)

    if(is.null(groups_mat)) {
      groups_mat <- sapply(unique(groups), function(x) x==groups)
      groups_mat <- cbind(groups_mat, -groups_mat)[, 1:ncol(groups_mat)]
    }

    # groups_mat <- cbind(-groups_mat, groups_mat)

  } else {

    groups_mat <- NULL
  }


  if (all(dim(sigma) == 1)) { return(1) }

  if (!is.null(mu)) {
    dvec <- mu
  } else {
    dvec <- rep(0, n)
  }

    # With return target
    safeOpt <- purrr::safely(quadprog::solve.QP)

    # Deprecated approach:
    # Using gamma as a hard constraint as below --> produces corner solutions.
    # This format is a Return-maximisation with a risk penalty, BUT with a cap on expected returns...
    # Amat <- cbind(1, -dvec, -diag(n), diag(n), groups_mat)
    # bvec <- c(1, -gamma, -UB, LB, -group.UB, group.LB)
    # opt <- safeOpt(sigma, dvec, Amat, bvec, meq = 1)

    # New approach:
    # The following fixes the objective tradeoff implicitly to be:
    # max μ′w − gamma x w′Σ w
    # gamma applied to covariance - so sensitivity is to cov, not return capping.
    # Gamma now a proper Risk aversion parameter, as opposed to a coinstraint on returns.
    # Small gamma → aggressive / return-seeking, big gamma: conservative
    # Thus: U(w) = μ′w − γ x σ x 2(w)

    Amat <- cbind( 1, -diag(n), diag(n), -groups_mat, groups_mat)
    bvec <- c( 1,-UB,LB,-group.UB,group.LB)

    if (gamma <= 0) stop("gamma must be > 0")
    Dmat <- 2 * sigma
    dvec_use <- dvec / gamma

    # This is gives us: max (mu′w − gamma w ′ Σ w)
    # QP minimises 0.5 x'Dmat x - dvec'x:

    opt <- safeOpt(Dmat      = Dmat,
                   dvec      = dvec_use,
                   Amat      = con$Amat,
                   bvec      = con$bvec,
                   meq       = 1)

    if(!is.null(opt$error)) return(stop(glue::glue("Optimization Error:\n\n {opt$error} \n\n")))
    opt_weights <- opt$result$solution


  if(!is.null(opt_weights)) {names(opt_weights) <- asset_names}
  if(!is.null(opt_weights) & !is.null(groups_mat)) {
    if(!all(opt_weights %*% groups_mat <= group.UB+1e-8)) stop("Violation of UB")
    if(!all(opt_weights %*% groups_mat >= group.LB-1e-8)) stop("Violation of LB")
  }

  return(opt_weights)

}

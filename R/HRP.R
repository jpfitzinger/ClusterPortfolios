#' @name HRP
#' @title Hierarchical Risk Parity portfolio
#' @description Computes optimal HRP portfolio with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix and does not need to be positive definite.
#'
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#'
#' By default, HRP constructs a balanced graph by bisecting each cluster at the
#' central location \eqn{(tau = 0)}{(tau = 0)}. Alternatively it is possible to split along
#' the actual cluster dendrogram \eqn{(tau = 1)}{(tau = 1)}, or to trade-off the objectives
#' by setting \eqn{0 < tau < 1}{0 < tau < 1}. In the latter case, \code{tau} is the proportion
#' of nodes on the left and right of the central split from which to select the
#' optimal splitting location (i.e. split at the optimal location that is at
#' most \eqn{tau\times x N/2}{tau x N/2} elements from the central split).
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param tau trade-off between naive (0) or cluster-based (1) tree-splitting (see Details).
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author Johann Pfitzinger
#' @references
#' Lopez de Prado, M. (2016).
#' Building Diversified Portfolios that Outperform Out-of-Sample.
#' \emph{SSRN Electronic Journal}.
#'
#' Pfitzinger, J., Katzke, N. (2019).
#' A Constrained Hierarchical Risk Parity Algorithm with Cluster-Based Capital Allocation.
#' \emph{Stellenbosch University, Department of Economics}. Working Paper 14/2019.
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' HRP(sigma, UB = 0.15, tau = 0.5)
#'
#' @export

HRP <- function(
  sigma,
  cluster_method = c("single", "average", "complete", "ward", "DIANA"),
  UB = NULL,
  LB = NULL,
  tau = 0
  ) {

  cluster_method <- match.arg(cluster_method)

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

  # Create Hierarchy
  S <- .compute_S_matrix(sigma, cluster_method = cluster_method, tau = tau)
  P <- .compute_P_matrix(S)

  w <- rep(NA, nrow(S))
  w.cumprod <- rep(NA, nrow(S))

  # Loop through nodes
  for (i in 1:nrow(S)) {

    if (is.na(w[i])) {

      # Identify adjacent nodes
      idx <- apply(P, 1, function(x) all(x == P[i,]))

      # Identify parent node
      parNode <- which(rowSumsC(as.matrix(S[, S[i,]==1])) == max(rowSumsC(as.matrix(S[, S[i,]==1]))))
      parNode <- parNode[rowSumsC(S)[parNode] > sum(S[i, ])]
      if (length(parNode)>0) parNode <- parNode[rowSumsC(S)[parNode] == min(rowSumsC(S)[parNode])]

      # Calculate node-specific constraints
      UBsub <- as.numeric(S %*% UB)[idx] / ifelse(length(parNode)==0, 1, w.cumprod[parNode])
      LBsub <- as.numeric(S %*% LB)[idx] / ifelse(length(parNode)==0, 1, w.cumprod[parNode])

      # Check if constraints are binding and set weights equal constraints to save an optimisation step
      if (round(sum(LBsub), 6) == 1) {

        w[idx] <- LBsub
        if (length(parNode) == 0) w.cumprod[idx] <- LBsub else w.cumprod[idx] <- w.cumprod[parNode] * LBsub

      } else if (round(sum(UBsub), 6) == 1) {

        w[idx] <- UBsub
        if (length(parNode) == 0) w.cumprod[idx] <- UBsub else w.cumprod[idx] <- w.cumprod[parNode] * UBsub

      } else {

        # Calculate cluster variance
        alpha <- sapply(which(idx), function(x) {
          ivp <- 1 / diag(as.matrix(sigma[S[x,] == 1, S[x,] == 1]))
          ivp <- ivp / sum(ivp)
          ivp <- ivp %*% sigma[S[x,] == 1, S[x,] == 1] %*% ivp
        })

        # Note: LdP does not take the sqrt() of alpha! This is for comparison with invvol
        alpha <- 1 / alpha / sum(1 / alpha)

        # Resolve constraint violations
        # The idea is simple: if constraint is violated, set to binding limit (i.e. make equal to constraint)
        # Subsequently distribute any excess allocation to the remaining weights in proportion of their inverse variance
        delta <- pmax(0, LBsub - alpha) - pmax(0, alpha - UBsub)
        maxit <- 100
        niter <- 0
        while (any(abs(delta) > 0) && niter < maxit) {

          alpha[abs(delta) > 0] <- alpha[abs(delta) > 0] + delta[abs(delta) > 0]
          alpha[delta == 0] <- alpha[delta == 0] + (1 - sum(alpha)) * alpha[delta == 0] / sum(alpha[delta == 0])
          alpha <- alpha / sum(alpha)
          delta <- pmax(0, LBsub - alpha) - pmax(0, alpha - UBsub)
          niter <- niter + 1

        }

        w[idx] <- alpha
        if (length(parNode) == 0) w.cumprod[idx] <- alpha else w.cumprod[idx] <- w.cumprod[parNode] * alpha

      }

    }

  }

  # Collapse weights to the asset dimension
  w.out <- c()
  for (i in 1:ncol(S)) w.out[i] <- prod(w[S[,i]==1])

  w.out[idx.short] <- -w.out[idx.short]

  names(w.out) <- asset_names

  return(w.out)
}

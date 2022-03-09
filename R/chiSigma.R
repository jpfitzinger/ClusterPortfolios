#' @name chiSigma
#' @title A convex hierarchical filter of the covariance matrix
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
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of estimated returns.
#' @param meta_loss a loss function of the most diversified hierarchical allocation graph.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param gamma risk aversion parameter. Default: \code{gamma = 0}.
#' @param max_tilt maximum percentage reduction in the effective number of assets. Default: \code{max_tilt = 1} (no restriction).
#' @param ... arguments passed to \code{cluster::agnes} method.
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
#'
#' @importFrom abind abind
#' @importFrom nloptr slsqp
#' @importFrom Matrix bdiag
#' @importFrom MASS ginv
#' @importFrom purrr quietly


chiSigma <- function(
  sigma,
  mu = NULL,
  meta_loss = c("MaxDiv", "ERC"),
  UB = NULL,
  LB = NULL,
  gamma = 0,
  max_tilt = 1,
  ...
) {

  meta_loss <- match.arg(meta_loss)
  quiet_slsqp <- purrr::quietly(nloptr::slsqp)

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

  cluster_object <- .get_clusters(sigma, ...)

  lvl_w <- function(n_clusters, sigma, cluster_object) {

    cut <- cutree(cluster_object, k = n_clusters)
    max_cut <- max(cut)

    cut_fx <- function(rowSel,cut) as.data.frame(matrix(as.numeric(rowSel == cut), ncol = length(cut)))
    S_Filler <- lapply(1:max_cut, cut_fx, cut = cut)
    S = matrix(nrow = length(S_Filler), ncol = length(cut))
    # Can be improved with Rcpp if required.
    for(i in 1:length(S_Filler) ) S[i,] <- as.matrix(S_Filler[i][[1]])

    S_av <- sweep(S, 1, rowSums(S), "/")
    sigma_sub <- S_av %*% sigma %*% t(S_av)

    if (!is.null(mu)) {
      mu_sub <- drop(mu %*% t(S_av))
    } else {
      mu_sub <- rep(0, max_cut)
    }

    UB_sub <- drop(UB %*% t(S))
    LB_sub <- drop(LB %*% t(S))

    Amat <- cbind(1, -diag(max_cut), diag(max_cut))
    bvec <- c(1, -UB_sub, LB_sub)

    opt <- quadprog::solve.QP(sigma_sub, mu_sub * gamma, Amat, bvec, meq = 1)
    w <- as.numeric(t(S_av) %*% opt$solution)

    sigma_ <- cbind(rbind(sigma, 1), c(rep(1, nrow(sigma)), 0))
    S_av_ <- as.matrix(Matrix::bdiag(S_av, 1))
    sigma_sub_ <- cbind(rbind(sigma_sub, 1), c(rep(1, nrow(sigma_sub)), 0))
    rot_mat <- t(S_av_) %*% solve(sigma_sub_) %*% S_av_ %*% sigma_

    return(list(w = w, rotation = rot_mat))

  }

  level_k <- c(n:1)
  w_opt_lvl <- lapply(level_k, lvl_w, sigma = sigma, cluster_object = cluster_object$cluster_object)
  w_mat <- sapply(w_opt_lvl, function(x) x$w)

  # Drop levels with same weights
  ix <- rep(T, ncol(w_mat))
  for (i in 2:ncol(w_mat)) {
    #ix[i] <- !isTRUE(all.equal(w_mat[,i-1], w_mat[,i], tolerance = 1e-8))
  }
  w_mat <- w_mat[, ix]
  n_meta <- ncol(w_mat)

  sigma_meta <- t(w_mat) %*% sigma %*% w_mat

  if (meta_loss == "MaxDiv") {

    .pRC <- function(w, w_mat, sigma_meta, sigma) {
      sigmaw <- crossprod(sigma_meta, w)
      w_ <- drop(w_mat %*% w)
      #pDR <- sqrt(as.numeric(crossprod(w, sigmaw))) / crossprod(w, sqrt(diag(sigma_meta)))
      pDR <- as.numeric(-crossprod(w_, sqrt(diag(sigma)))) / sqrt(as.numeric(crossprod(w, sigmaw)))
      return(pDR)
    }
    .eqConstraint <- function (w)
    {
      return(sum(w) - 1)
    }
    .hinConstraint <- function(w)
    {
      return(crossprod(level_k[ix], w) - ((n-1)*(1-max_tilt)+1))
    }
    phi <- quiet_slsqp(x0 = rep(1/n_meta, n_meta), fn = .pRC,
                         heq = .eqConstraint,
                         hin = .hinConstraint,
                         lower = rep(1e-5, n_meta),
                         upper = rep(1, n_meta), nl.info = FALSE,
                         control = list(xtol_rel = 1e-18, check_derivatives = FALSE,
                                        maxeval = 20000),
                         sigma_meta = sigma_meta, w_mat = w_mat, sigma = sigma)$result$par

  }

  if (meta_loss == "ERC") {

    .pRC <- function(w, w_mat, sigma) {
      w <- as.numeric(w_mat %*% w)
      sigmaw <- crossprod(sigma, w)
      pRC <- (w * sigmaw)/as.numeric(crossprod(w, sigmaw))
      d <- sum((pRC - 1/n_meta)^2)
      return(d)
    }
    .eqConstraint <- function (w)
    {
      return(sum(w) - 1)
    }
    .hinConstraint <- function(w)
    {
      return(crossprod(level_k[ix], w) - ((n-1)*(1-max_tilt)+1))
    }
    phi <- quiet_slsqp(x0 = rep(1/n_meta, n_meta), fn = .pRC,
                         heq = .eqConstraint,
                         hin = .hinConstraint,
                         lower = rep(1e-5, n_meta),
                         upper = rep(1, n_meta), nl.info = FALSE,
                         control = list(xtol_rel = 1e-18, check_derivatives = FALSE,
                                        maxeval = 20000),
                         sigma = sigma, w_mat = w_mat)$result$par

  }

  rot_mat <- lapply(w_opt_lvl, function(x) x$rotation)
  rot_mat <- abind::abind(rot_mat[ix], along = 3)
  rot_mat <- apply(rot_mat, 2, function(x) x %*% phi)
  rot_mat <- MASS::ginv(rot_mat[1:n, 1:n])

  chiSigma_Mat <- t(rot_mat) %*% sigma %*% rot_mat
  diag(chiSigma_Mat) <- diag(chiSigma_Mat) + 1e-8
  colnames(chiSigma_Mat) <- rownames(chiSigma_Mat) <- asset_names

  out <- list(sigma = chiSigma_Mat)

  if (!is.null(mu)) {
    chiMu_Vec <- drop(mu %*% rot_mat)
    names(chiMu_Vec) <- asset_names
    out$mu <- chiMu_Vec
  }

  out$phi <- phi
  out$w <- drop(w_mat %*% phi)

  return(out)

}


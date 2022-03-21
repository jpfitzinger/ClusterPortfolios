#' @name chiPlot
#' @title Plot the dendrogram of a CHI portfolio
#' @description Plots the dendrogram of a CHI portfolio.
#'
#' @details The dendrogram is generated using hierarchical clustering and modified
#' so that the height differential between any two splits is the shrinkage weight of
#' the lower split (ranging between \code{0} and \code{1}). With no shrinkage, all shrinkage weights
#' are equal to \code{1} and the dendrogram has a height of \eqn{p}{p}. With shrinkage
#' the dendrogram has a height of \eqn{(\kappa \times p)}{(\code{kappa} x p)}.
#'
#' The leaf nodes are colored to indicate the coefficient sign, with the size indicating
#' the absolute magnitude of the coefficients.
#'
#' A color bar on the right indicates the relative contribution of each level to the
#' coefficient of determination, with darker hues representing a larger contribution.
#'
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
#' @param max_leaf_size maximum size of the leaf nodes. Default is \code{max_leaf_size=3}.
#' @param ymax upper limit for y axis
#' @param ... arguments passed to \code{cluster::agnes} method.
#' @return A plotted dendrogram.
#' @author Johann Pfitzinger
#'
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' chiPlot(sigma, UB = 0.15)
#'
#' @export
#'
#' @seealso \code{CHI} and \code{chiSigma} methods

chiPlot <- function(
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
  max_leaf_size = 3,
  ymax = NULL,
  ...
) {

  meta_loss <- match.arg(meta_loss)

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  chi <- chiSigma(sigma = sigma, mu = mu, meta_loss = meta_loss,
                  UB = UB, LB = LB, gamma = gamma, max_tilt = max_tilt,
                  ...)

  w <- MV(sigma = chi$sigma, mu = chi$mu, UB = UB, LB = LB, gamma = gamma,
          groups = groups, group.UB = group.UB, group.LB = group.LB)

  clust <- chi$cluster_object
  phi <- chi$phi

  aggr <- diag(length(phi))
  aggr[lower.tri(aggr)] <- 1
  theta <- rev(as.numeric(aggr %*% phi))

  dof <- rev(chi$n_assets_per_level)
  dof <- c(dof[1], dof[-1] - dof[-length(dof)])

  heights <- theta * dof

  expl_variance <- rev(chi$expl_var)
  #expl_variance <- c(diff(expl_variance), 0)


  .draw_dendro(clust, w, heights, expl_variance, asset_names, x$df,
               max_leaf_size, ymax)

}

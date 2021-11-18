#' @name hSigma
#' @title Hierarchical filtering of the covariance matrix
#' @description Generates a hierarchically filtered covariance matrix than can be used for optimization.
#' @details The argument \code{sigma} is a covariance matrix.
#'
#' Hierarchical clustering is performed using the \code{cluster}-package. If
#' \code{cluster_method == 'DIANA'}, the function \code{cluster::diana} is used
#' to compute a cluster dendrogram, otherwise the function \code{cluster::agnes(., method = cluster_method)}
#' is used. Default is single-linkage agglomerative nesting.
#'
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param cluster_method hierarchical cluster algorithm used to construct an asset hierarchy.
#' @return A \eqn{(N \times N)}{(N x N)} filtered covariance matrix.
#' @author Johann Pfitzinger
#' @references
#'
#' @examples
#' # Load returns of assets or portfolios
#' data("Industry_10")
#' rets <- Industry_10
#' sigma <- cov(rets)
#' hsig <- hSigma(sigma)
#' MV(hsig)
#'
#' @export


hSigma <- function(
  sigma,
  cluster_method = c("single", "average", "complete", "ward", "DIANA")
  ) {

  cluster_method <- match.arg(cluster_method)

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  # Cluster
  clust <- .get_clusters(sigma, cluster_method, 2)
  cut <- clust$clusters
  clust <- clust$cluster_object

  # Create S list
  S <- purrr::map(1:n, function(k) {
    cut <- cutree(clust, k)
    max_cut <- max(cut)
    cut_fx <- function(rowSel,cut) as.data.frame(matrix(as.numeric(rowSel == cut), ncol = length(cut)))
    S_Filler <- purrr::map(1:max_cut, ~cut_fx(., cut))
    S = matrix(nrow = length(S_Filler), ncol = length(cut))
    # Can be improved with Rcpp if required.
    for(i in 1:length(S_Filler) ) S[i,] <- as.matrix(S_Filler[i][[1]])
      return(S)
    })

  S_av <- lapply(S, function(S) {
    S_av <- sweep(S, 1, rowSums(S), "/")
  })

  # For each S matrix, filter
  filtered_corr <- matrix(1, n, n)
  for (i in 1:n) {

    fil_inner <- t(S[[i]]) %*% S_av[[i]] %*% cov2cor(sigma) %*% t(S_av[[i]]) %*% S[[i]]
    ix <- round(cov2cor(filtered_corr), 6)==1
    filtered_corr[ix] <- fil_inner[ix]

  }

  filtered_cov <- diag(sqrt(diag(sigma))) %*% filtered_corr %*% diag(sqrt(diag(sigma)))
  colnames(filtered_cov) <- rownames(filtered_cov) <- asset_names

  return(filtered_cov)

}


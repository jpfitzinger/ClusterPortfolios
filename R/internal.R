#' @importFrom cluster agnes
#' @importFrom cluster diana
#' @importFrom cluster clusGap
#' @importFrom stats cov2cor
#' @importFrom stats dist
#' @importFrom stats cutree

.cluster_object <- function(sigma, cluster_method) {

  corr <- cov2cor(sigma)
  #distmat <- ((1 - corr) / 2)^0.5
  distmat <- dist(corr)
  if (cluster_method == "DIANA") {
    clust <- diana(as.dist(distmat))
  } else {
    clust <- agnes(as.dist(distmat), method = cluster_method)
  }

  return(clust)

}

.get_clusters <- function(sigma, cluster_method, n_clusters = NULL) {

  clust <- .cluster_object(sigma, cluster_method)

  if (is.null(n_clusters)) {

    clusters <- NULL

  } else if (n_clusters == "auto") {

    if (cluster_method == "DIANA") {
      warning("'DIANA' clustering currently not supported for automatic cluster selection. Changing cluster_method to 'complete'.")
      cluster_method <- "complete"
    }
    max_clusters <- ceiling(dim(sigma)[1] / 2)
    corr <- cov2cor(sigma)
    distmat <- ((1 - corr) / 2)^0.5
    opt_clust_fx <- purrr::quietly(NbClust::NbClust)
    opt <- opt_clust_fx(diss = as.dist(distmat), max.nc = max_clusters,
                        method = cluster_method, index = "silhouette", distance = NULL)$result
    opt_n_clusters <- opt$Best.nc[1]
    clusters <- factoextra::hcut(as.dist(distmat), k = opt_n_clusters,
                                 hc_func = "agnes",
                                 hc_method = cluster_method, isdiss = T)$cluster

  } else if (is.integer(n_clusters) | is.numeric(n_clusters)) {

    clusters <- cutree(clust, k = n_clusters)

  }

  return(list(cluster_object = clust, clusters = clusters))

}

.compute_S_matrix <- function(sigma, cluster_method, tau) {

  clust <- .get_clusters(sigma, cluster_method)

  cluster_order <- clust$cluster_object$order
  cluster_height <- clust$cluster_object$height

  bisect.inner <- function(cluster_order, cluster_height, tau, corr, sigma) {

    N <- length(cluster_order)
    idx <- floor(N/2 - (N/2 - 1) * tau):floor(N/2 + (N/2 - 1) * tau)
    split.pos <- idx[which.max(cluster_height[idx])]

    x.vec.1 <- 1:split.pos
    x.vec.2 <- c(1:N)[-x.vec.1]

    vec.1 <- vec.2 <- rep(0, ncol(sigma))
    vec.1[cluster_order[x.vec.1]] <- 1
    vec.2[cluster_order[x.vec.2]] <- 1

    S.mat <<- cbind(S.mat, vec.1, vec.2)

    if (length(x.vec.1) > 1) bisect.inner(cluster_order[x.vec.1], cluster_height[x.vec.1], tau, corr[x.vec.1, x.vec.1], sigma)
    if (length(x.vec.2) > 1) bisect.inner(cluster_order[x.vec.2], cluster_height[x.vec.2], tau, corr[x.vec.2, x.vec.2], sigma)

  }

  S.mat <- c()

  bisect.inner(cluster_order, cluster_height, tau, corr, sigma)

  return(t(S.mat))

}

.compute_P_matrix <- function(S) {

  Sfull <- rbind(1, S)

  SparentSums <- sapply(1:nrow(S), function(x) rowSumsC(as.matrix(Sfull[, S[x,]==1])))
  Sparent <- S

  for (i in 1:nrow(S)) {

    idx <- SparentSums[i + 1, i] == SparentSums[, i]
    idx[i + 1] <- F

    Ssize <- rowSumsC(Sfull)
    Ssize[!idx] <- max(Ssize) + 1

    Sparent[i,] <- Sfull[which.min(Ssize),]

  }

  P <- PMat(Sparent)

  return(P)

}

#' @importFrom cluster agnes
#' @importFrom cluster diana
#' @importFrom stats cov2cor
#' @importFrom stats dist

.compute_S_matrix <- function(sigma, cluster_method, tau) {

  corr <- cov2cor(sigma)
  distmat <- ((1 - corr) / 2)^0.5
  if (cluster_method == "DIANA") {
    clust <- cluster::diana(dist(distmat))
  } else {
    clust <- cluster::agnes(dist(distmat), method = cluster_method)
  }

  cluster_order <- clust$order
  cluster_height <- clust$height

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

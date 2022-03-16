#' @importFrom cluster agnes silhouette
#' @importFrom stats cov2cor dist cutree as.dist

.cluster_object <- function(sigma, ...) {

  corr <- stats::cov2cor(sigma)
  distmat <- stats::dist(corr)
  clust <- cluster::agnes(stats::as.dist(distmat), ...)

  return(clust)

}

.get_clusters <- function(sigma, n_clusters = NULL, ...) {

  clust <- .cluster_object(sigma, ...)

  if (is.null(n_clusters)) {

    clusters <- NULL

  } else if (n_clusters == "auto") {

    max_clusters <- ceiling(dim(sigma)[1] / 2)
    corr <- stats::cov2cor(sigma)
    distmat <- stats::dist(corr)
    clust <- .cluster_object(sigma, ...)

    cluster_no <- c(2:max_clusters)

    sil_widths <- sapply(cluster_no, function(k, clust, distmat) {
      mean(cluster::silhouette(stats::cutree(clust, k = k), dist = distmat)[, 3])
    }, clust = clust, dist = distmat)

    opt_n_clusters <- cluster_no[which.max(sil_widths)]
    clusters <- stats::cutree(clust, k = opt_n_clusters)

  } else if (is.integer(n_clusters) | is.numeric(n_clusters)) {

    clusters <- cutree(clust, k = n_clusters)

  }

  return(list(cluster_object = clust, clusters = clusters))

}

.compute_S_matrix <- function(sigma, tau, ...) {

  clust <- .get_clusters(sigma, ...)

  cluster_order <- clust$cluster_object$order
  cluster_height <- clust$cluster_object$height

  bisect.inner <- function(cluster_order, cluster_height, tau, sigma) {

    N <- length(cluster_order)
    idx <- floor(N/2 - (N/2 - 1) * tau):floor(N/2 + (N/2 - 1) * tau)
    split.pos <- idx[which.max(cluster_height[idx])]

    x.vec.1 <- 1:split.pos
    x.vec.2 <- c(1:N)[-x.vec.1]

    vec.1 <- vec.2 <- rep(0, ncol(sigma))
    vec.1[cluster_order[x.vec.1]] <- 1
    vec.2[cluster_order[x.vec.2]] <- 1

    S.mat <<- cbind(S.mat, vec.1, vec.2)

    if (length(x.vec.1) > 1) bisect.inner(cluster_order[x.vec.1], cluster_height[x.vec.1], tau, sigma)
    if (length(x.vec.2) > 1) bisect.inner(cluster_order[x.vec.2], cluster_height[x.vec.2], tau, sigma)

  }

  S.mat <- c()

  bisect.inner(cluster_order, cluster_height, tau, sigma)

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

.draw_dendro <- function(clust, w, heights, explained_variance,
                         asset_names, df, max_leaf_size) {

  w_sizes <- abs(w)/max(abs(w))
  w_sizes <- w_sizes * max_leaf_size
  w_col <- ifelse(sign(w)>=0, "darkblue", "firebrick")

  n <- length(heights)
  dend_heights <- cumsum(rev(heights))
  clust$height <- dend_heights[-n]

  dend <- stats::as.dendrogram(clust)
  dend <- dendextend::set(dend, "labels", asset_names[clust$order])
  dend <- dendextend::set(dend, "leaves_pch", 15)
  dend <- dendextend::set(dend, "leaves_cex", w_sizes[clust$order])
  dend <- dendextend::set(dend, "leaves_col", w_col[clust$order])

  pal <- RColorBrewer::brewer.pal(9, "Blues")
  pal <- c("#FFFFFF", pal)
  cols <- rev(explained_variance)
  cols <- pmax(c(cols[-n] - cols[-1], cols[n]), 0)
  cols <- round(sqrt((cols - min(cols)) / (max(cols) - min(cols))) * (length(pal)-1)+1)

  top_node <- dendextend::get_nodes_xy(dend)[1,]

  graphics::plot(stats::as.dendrogram(dend), ylim=c(0, max(dend_heights)),
                 panel.first=graphics::abline(h = dend_heights[dend_heights > 1e-4],
                                              col="lightgrey", lwd=1, lty = "dashed"))
  graphics::segments(x0 = top_node[1], y0 = top_node[2], y1 = dend_heights[n])
  graphics::points(x = top_node[1], y = dend_heights[n], pch = 15)
  #graphics::rect(n*1.015, c(0, dend_heights[-n]), n*1.03, dend_heights, col = pal[cols], lwd=0.1)

}

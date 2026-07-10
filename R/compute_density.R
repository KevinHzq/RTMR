#' Kernel density of point features at grid cell centroids
#'
#' Computes an Epanechnikov kernel density estimate of the point features,
#' evaluated at each grid cell centroid, as done by the RTMDx Utility for the
#' "density" operationalization of spatial influence.
#'
#' @param pt An sf or sfc object with the point features of a risk factor.
#' @param grid An sf grid created by [create_grid()].
#' @param bandwidth Kernel bandwidth (search radius) in the units of the CRS.
#'
#' @return A numeric vector of density values, one per grid cell.
#' @export
compute_density <- function(pt, grid, bandwidth) {
  grid_centroid <- suppressWarnings(sf::st_centroid(sf::st_geometry(grid)))
  pt_geom <- sf::st_geometry(pt)

  # only point/centroid pairs within the bandwidth contribute to the kernel sum
  hits <- sf::st_is_within_distance(grid_centroid, pt_geom, dist = bandwidth)

  density <- numeric(length(grid_centroid))

  cell_idx <- rep(seq_along(hits), lengths(hits))
  if (length(cell_idx) == 0) {
    return(density)
  }
  pt_idx <- unlist(hits)

  d <- as.numeric(
    sf::st_distance(grid_centroid[cell_idx], pt_geom[pt_idx], by_element = TRUE)
  )

  # 2D Epanechnikov kernel: K(u) = 2 / (pi * h^2) * (1 - u^2), u = d / h <= 1
  k <- 2 / (pi * bandwidth^2) * (1 - (d / bandwidth)^2)

  contrib <- rowsum(k, group = cell_idx)
  density[as.integer(rownames(contrib))] <- contrib[, 1]

  density
}

#' Binarize density values into high-density cells
#'
#' Reclassifies kernel density values into high-density (1) and other (0)
#' cells. Following RTMDx, high-density cells are by default those with
#' values at least `n_sd` (2) standard deviations above the mean density.
#' A constant density surface (e.g. no features within reach of any cell)
#' has no high-density cells, so all zeros are returned.
#'
#' The mean + 2 SD rule flags about the top 2% of cells when densities are
#' roughly normal, but kernel densities of point patterns are strongly
#' right-skewed, so the flagged share drifts uncontrollably with the skew
#' — and hence differs between risk factors and bandwidths. Supplying
#' `quantile` switches to a quantile rule that pins the flagged share
#' directly: `quantile = 0.95` flags the cells above the 95th percentile
#' (the top 5%), making the "high-density" exposure comparable across
#' factors. The threshold is applied strictly (`>`), so the flagged share
#' is at most `1 - quantile`; ties can only reduce it. In particular, when
#' more than `quantile` of the cells have zero density (out of reach of
#' every feature), the threshold lands on 0 and exactly the
#' positive-density cells are flagged, and a constant surface flags
#' nothing.
#'
#' @param x Numeric vector of density values, e.g. from [compute_density()].
#' @param n_sd Number of standard deviations above the mean defining the
#'   high-density threshold (default 2, as in RTMDx). Ignored when
#'   `quantile` is given.
#' @param quantile Optional quantile cutoff strictly between 0 and 1;
#'   cells with density above this quantile are high-density. 0.95 is a
#'   sensible choice (top 5% of cells).
#'
#' @return An integer vector of 0/1 values.
#' @export
binarize_density <- function(x, n_sd = 2, quantile = NULL) {
  if (!is.null(quantile)) {
    if (!is.numeric(quantile) || length(quantile) != 1 ||
      quantile <= 0 || quantile >= 1) {
      stop("`quantile` must be a single number strictly between 0 and 1")
    }
    thr <- stats::quantile(x, quantile, names = FALSE)
    return(as.integer(x > thr))
  }
  s <- stats::sd(x)
  if (is.na(s) || s == 0) {
    return(integer(length(x)))
  }
  as.integer(x >= mean(x) + n_sd * s)
}

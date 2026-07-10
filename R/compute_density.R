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
#' cells. Following RTMDx, high-density cells are those with values at least
#' 2 standard deviations above the mean density.
#'
#' @param x Numeric vector of density values, e.g. from [compute_density()].
#' @param n_sd Number of standard deviations above the mean defining the
#'   high-density threshold (default 2).
#'
#' @return An integer vector of 0/1 values.
#' @export
binarize_density <- function(x, n_sd = 2) {
  as.integer(x >= mean(x) + n_sd * stats::sd(x))
}

#' Distance from each grid cell to the nearest feature
#'
#' Computes the Euclidean distance from each grid cell centroid to the
#' nearest feature. This is the basis of the "proximity" operationalization
#' of a risk factor's spatial influence. Features may be points, lines, or
#' polygons: for lines and polygons the distance is to the nearest edge,
#' and 0 if the centroid falls inside a polygon.
#'
#' @param pt An sf or sfc object with the features of a risk factor.
#' @param grid An sf grid created by [create_grid()].
#'
#' @return A numeric vector of distances (in the units of the CRS), one per
#'   grid cell.
#' @export
compute_proximity <- function(pt, grid) {
  grid_centroid <- suppressWarnings(sf::st_centroid(sf::st_geometry(grid)))
  pt_geom <- sf::st_geometry(pt)

  nearest_pt <- sf::st_nearest_feature(grid_centroid, pt_geom)

  as.numeric(sf::st_distance(grid_centroid, pt_geom[nearest_pt], by_element = TRUE))
}

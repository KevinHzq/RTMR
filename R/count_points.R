# count number of points in each grid
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r

#' Count points falling in each grid cell
#'
#' Typically used to build the outcome variable (e.g. crime counts per cell).
#'
#' @param pt An sf or sfc object with point features (e.g. crime events).
#' @param grid An sf grid created by [create_grid()].
#'
#' @return An integer vector of counts, one per grid cell.
#' @export
count_points <- function(pt, grid) {
  lengths(sf::st_intersects(grid, pt))
}

#' Count points falling in each grid cell
#'
#' Typically used to build the outcome variable (e.g. crime counts per cell).
#'
#' Each point is counted exactly once: a point lying on an edge shared by
#' two cells is assigned to the first intersecting cell in grid order (a
#' fixed half-open convention, as in raster software). Points outside the
#' grid are dropped.
#'
#' @param pt An sf or sfc object with point features (e.g. crime events).
#' @param grid An sf grid created by [create_grid()].
#'
#' @return An integer vector of counts, one per grid cell.
#' @export
count_points <- function(pt, grid) {
  hits <- sf::st_intersects(sf::st_geometry(pt), grid)
  cell <- vapply(
    hits,
    function(i) if (length(i) > 0) i[1] else NA_integer_,
    integer(1)
  )
  tabulate(cell[!is.na(cell)], nbins = length(sf::st_geometry(grid)))
}

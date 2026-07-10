#' Create an analysis grid over a study area
#'
#' Builds a regular grid of cells (the units of analysis) covering the study
#' area. Cells are kept if they have at least some overlap with the study
#' area polygon, mirroring the masking behaviour of the RTMDx Utility.
#'
#' @param x An sf object with the study area polygon(s).
#' @param ... Passed on to [sf::st_make_grid()], e.g. `cellsize`,
#'   `square = FALSE` for hexagons.
#' @param clip If `TRUE` (default), keep only cells intersecting the study
#'   area. If `FALSE`, return the full bounding-box grid.
#'
#' @return An sf object with a `cell_id` column and polygon geometry.
#' @export
#'
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' grid <- create_grid(nc, n = c(20, 10))
create_grid <- function(x, ..., clip = TRUE) {
  xboundary <- sf::st_union(sf::st_geometry(x))

  sfc_grid <- sf::st_make_grid(xboundary, ...)

  sf_grid <- sf::st_sf(cell_id = seq_along(sfc_grid), geometry = sfc_grid)

  if (clip) {
    # equivalent to sf::st_filter(), which is avoided because it requires
    # dplyr (an optional dependency of sf)
    sf_grid[lengths(sf::st_intersects(sf_grid, xboundary)) > 0, ]
  } else {
    sf_grid
  }
}

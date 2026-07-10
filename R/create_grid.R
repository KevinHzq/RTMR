#' Create an analysis grid over a study area
#'
#' Builds a regular grid of cells (the units of analysis) covering the study
#' area. Cells are kept if they have at least some overlap with the study
#' area polygon, mirroring the masking behaviour of the RTMDx Utility.
#'
#' Areas where the outcome cannot occur (water bodies, restricted land)
#' can be removed with `exclude`: cells lying *entirely* within the
#' exclusion polygons are dropped. Such cells are structural zeros — they
#' would otherwise dilute the model with places that carry no information
#' about risk — and removing them restricts the analysis (and the risk
#' map) to places where events are possible. Cells only partially covered
#' by the exclusion are kept, since events can still occur in their
#' uncovered part.
#'
#' @param x An sf object with the study area polygon(s).
#' @param ... Passed on to [sf::st_make_grid()], e.g. `cellsize`,
#'   `square = FALSE` for hexagons.
#' @param clip If `TRUE` (default), keep only cells intersecting the study
#'   area. If `FALSE`, return the full bounding-box grid.
#' @param exclude Optional sf or sfc polygon layer of areas to exclude;
#'   cells fully covered by it are removed.
#'
#' @return An sf object with a `cell_id` column and polygon geometry.
#' @export
#'
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' grid <- create_grid(nc, n = c(20, 10))
create_grid <- function(x, ..., clip = TRUE, exclude = NULL) {
  xboundary <- sf::st_union(sf::st_geometry(x))

  sfc_grid <- sf::st_make_grid(xboundary, ...)

  sf_grid <- sf::st_sf(cell_id = seq_along(sfc_grid), geometry = sfc_grid)

  if (clip) {
    # equivalent to sf::st_filter(), which is avoided because it requires
    # dplyr (an optional dependency of sf)
    sf_grid <- sf_grid[lengths(sf::st_intersects(sf_grid, xboundary)) > 0, ]
  }

  if (!is.null(exclude)) {
    if (!inherits(exclude, c("sf", "sfc"))) {
      stop("`exclude` must be an sf or sfc object")
    }
    excl <- sf::st_union(sf::st_geometry(exclude))
    covered <- lengths(sf::st_covered_by(sf_grid, excl)) > 0
    if (all(covered)) {
      stop("`exclude` covers every grid cell; nothing is left to model")
    }
    sf_grid <- sf_grid[!covered, ]
  }

  sf_grid
}

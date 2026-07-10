#' Convert a risk terrain model to a raster
#'
#' Rasterizes the grid of a fitted [rtm()] model so the risk terrain map can
#' be used with raster tools (map algebra, comparison of risk surfaces) or
#' exported to GIS software. Cells outside the study area are `NA`, matching
#' the masking behaviour of the RTMDx Utility's GeoTiff outputs.
#'
#' Requires the terra package. The raster resolution is derived from the
#' grid cell area, which is exact for the default square grid; hexagonal
#' grids are rasterized at an approximately equivalent resolution.
#'
#' @param x An `rtm` object.
#' @param what Cell value to rasterize: `"relrisk"` (relative risk score,
#'   default; RTMDx's "output_score"), `"prediction"` (expected event count;
#'   RTMDx's "output_prediction"), or `"outcome_count"`.
#'
#' @return A [terra::SpatRaster] with one layer.
#' @seealso [write_rtm()] to write directly to a GeoTiff file.
#' @export
rtm_raster <- function(x, what = c("relrisk", "prediction", "outcome_count")) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("rtm_raster() requires the 'terra' package; install it with install.packages(\"terra\")")
  }
  what <- match.arg(what)

  cell_size <- sqrt(mean(as.numeric(sf::st_area(x$grid))))

  v <- terra::vect(x$grid[what])
  r <- terra::rast(terra::ext(v), resolution = cell_size, crs = terra::crs(v))
  terra::rasterize(v, r, field = what)
}

#' Write a risk terrain map to a GeoTiff
#'
#' Exports the risk terrain map of a fitted [rtm()] model as a GeoTiff (or
#' any other format supported by [terra::writeRaster()]), for use in GIS
#' software such as ArcGIS or QGIS. This mirrors the "output_score" and
#' "output_prediction" GeoTiffs produced by the RTMDx Utility.
#'
#' @param x An `rtm` object.
#' @param filename Output file path, e.g. `"risk_map.tif"`.
#' @param what Cell value to export; see [rtm_raster()].
#' @param ... Passed on to [terra::writeRaster()], e.g. `overwrite = TRUE`.
#'
#' @return The written [terra::SpatRaster], invisibly.
#' @export
write_rtm <- function(x, filename,
                      what = c("relrisk", "prediction", "outcome_count"),
                      ...) {
  r <- rtm_raster(x, what)
  terra::writeRaster(r, filename, ...)
  invisible(r)
}

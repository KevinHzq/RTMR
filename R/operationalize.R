#' Spatial influence distances to test for a risk factor
#'
#' Builds the sequence of threshold distances (for proximity) or kernel
#' bandwidths (for density) at which a risk factor's spatial influence is
#' tested, in half- or whole-block increments up to a maximum number of
#' blocks, as in the RTMDx Utility.
#'
#' @param block_length Average block length of the study area, in the units
#'   of the CRS (e.g. 500 for 500 ft, or 120 for 120 m).
#' @param max_blocks Maximum spatial influence, in blocks (RTMDx allows up to
#'   4; default 3).
#' @param increment `"whole"` (default) tests whole-block distances only;
#'   `"half"` also tests half-block distances.
#'
#' @return A numeric vector of distances.
#' @export
#'
#' @examples
#' spatial_influences(500, max_blocks = 3, increment = "half")
#' # 250 500 750 1000 1250 1500
spatial_influences <- function(block_length,
                               max_blocks = 3,
                               increment = c("whole", "half")) {
  increment <- match.arg(increment)
  step <- if (increment == "half") block_length / 2 else block_length
  seq(step, block_length * max_blocks, by = step)
}

#' Operationalize a risk factor into binary spatial-influence variables
#'
#' For one risk factor (a point layer), builds the set of candidate binary
#' variables tested by RTM: for each spatial influence distance, a
#' "proximity" variable (1 if the cell centroid is within the distance of a
#' feature, 0 otherwise) and/or a "density" variable (1 if the cell is in a
#' high-density area, i.e. the Epanechnikov kernel density with that
#' bandwidth is at least 2 standard deviations above the mean).
#'
#' Protective factors are coded 0/-1 instead of 0/1 (as in RTMDx) so that
#' model coefficients are expected to be positive for all variables.
#'
#' Line and polygon features are supported for the proximity
#' operationalization: distance is measured from the cell centroid to the
#' nearest edge of the feature, and cells whose centroid falls inside a
#' polygon are at distance 0 (always exposed). The density
#' operationalization is a kernel density of a *point pattern* and requires
#' point features; requesting it for other geometries is an error rather
#' than a misleading result (convert such features to representative points
#' first if a density reading is genuinely wanted).
#'
#' @param pt An sf or sfc object with the features of the risk factor
#'   (points for any operation; lines or polygons for proximity only).
#'   Include features slightly beyond the study area boundary so that their
#'   spatial influence on edge cells is captured.
#' @param grid An sf grid created by [create_grid()].
#' @param name Name of the risk factor, used to label the variables.
#' @param block_length,max_blocks,increment Passed to [spatial_influences()].
#' @param operation `"proximity"`, `"density"`, or `"both"` (default
#'   `"proximity"`, as in RTMDx).
#' @param type `"aggravating"` (default) if the factor is expected to
#'   increase risk, `"protective"` if expected to reduce it.
#'
#' @return A data.frame with one 0/1 (or 0/-1) column per candidate variable
#'   and one row per grid cell. The variable metadata (factor name,
#'   operation, spatial influence, type) is attached as attribute
#'   `"rtm_meta"`, a data.frame with one row per variable.
#' @export
operationalize <- function(pt, grid, name,
                           block_length,
                           max_blocks = 3,
                           increment = c("whole", "half"),
                           operation = c("proximity", "density", "both"),
                           type = c("aggravating", "protective")) {
  operation <- match.arg(operation)
  increment <- match.arg(increment)
  type <- match.arg(type)

  if (length(sf::st_geometry(pt)) == 0) {
    stop(
      "risk factor `", name, "` has no features; ",
      "check the layer and that it shares the grid's CRS"
    )
  }

  is_point <- all(as.character(sf::st_geometry_type(pt)) %in% c("POINT", "MULTIPOINT"))
  if (!is_point && operation %in% c("density", "both")) {
    stop(
      "the density operationalization requires point features; ",
      "use operation = \"proximity\" for line/polygon factors ",
      "or convert the features to representative points first"
    )
  }

  name <- make.names(name)
  distances <- spatial_influences(block_length, max_blocks, increment)
  sign_val <- if (type == "protective") -1L else 1L

  vars <- list()
  meta <- list()

  add_var <- function(op, d, values) {
    v <- paste0(name, "_", substr(op, 1, 4), "_", format(d, trim = TRUE, scientific = FALSE))
    vars[[v]] <<- sign_val * values
    meta[[v]] <<- data.frame(
      variable = v, factor = name, operation = op,
      spatial_influence = d, type = type,
      stringsAsFactors = FALSE
    )
  }

  if (operation %in% c("proximity", "both")) {
    dist_nearest <- compute_proximity(pt, grid)
    for (d in distances) {
      add_var("proximity", d, as.integer(dist_nearest <= d))
    }
  }

  if (operation %in% c("density", "both")) {
    for (d in distances) {
      add_var("density", d, binarize_density(compute_density(pt, grid, bandwidth = d)))
    }
  }

  out <- as.data.frame(vars, optional = TRUE)
  attr(out, "rtm_meta") <- do.call(rbind, c(meta, make.row.names = FALSE))
  out
}

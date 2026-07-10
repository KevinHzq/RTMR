#' Smooth denominator events into a total-preserving offset surface
#'
#' Builds a smoothed per-cell denominator for [rtm()] rate models by kernel
#' density estimation of the denominator events ([compute_density()],
#' Epanechnikov kernel), rescaled so that the cell values sum to the number
#' of denominator events falling in the grid (total-preserving). The result
#' can be passed to the `offset` argument of [rtm()] in place of raw
#' per-cell counts.
#'
#' Intended for *sensitivity analysis*, not as the primary offset. Raw
#' counts are statistically exact: cells with a zero denominator carry no
#' information about the event rate and their exclusion is harmless. A
#' smoothed offset instead (1) assigns synthetic denominators to cells with
#' no observed activity, (2) makes results sensitive to the `bandwidth`
#' choice, and (3) blurs denominator mass across exposure boundaries while
#' the outcome counts stay point-accurate, which can bias the RRVs of
#' factors co-located with activity peaks. Recommended use: fit the model
#' with raw counts as the primary analysis, refit with smoothed offsets at
#' two or three bandwidths, and report the stability of the RRVs.
#'
#' @param pt An sf or sfc object with the denominator event points (e.g.
#'   all overdose events).
#' @param grid An sf grid created by [create_grid()].
#' @param bandwidth Kernel bandwidth in the units of the CRS.
#'
#' @return A numeric vector with one non-negative value per grid cell,
#'   summing to the number of denominator events within the grid.
#' @seealso [rtm()], [compute_density()]
#' @export
smooth_offset <- function(pt, grid, bandwidth) {
  dens <- compute_density(pt, grid, bandwidth)
  total <- sum(count_points(pt, grid))

  if (total == 0 || sum(dens) == 0) {
    stop("no denominator events fall within reach of the grid; check the CRS and bandwidth")
  }

  dens * total / sum(dens)
}

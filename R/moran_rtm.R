#' Test model residuals for spatial autocorrelation with Moran's I
#'
#' Computes Moran's I on the residuals of a fitted [rtm()] model and tests
#' it against the null hypothesis of no spatial autocorrelation by random
#' permutation of the residuals across cells.
#'
#' The RTMDx procedure (and hence [rtm()]) treats grid cells as independent
#' observations. Heffner (2013) acknowledges that spatial autocorrelation is
#' not modeled and that significance values are inaccurate in its presence;
#' left unchecked it also makes the cross-validated culling optimistic and
#' weakens BIC's complexity penalty, both of which push towards selecting
#' too many variables. This diagnostic makes the assumption checkable:
#' residual spatial clustering that the selected risk factors have *not*
#' absorbed shows up as a positive Moran's I with a small permutation
#' p-value.
#'
#' Neighbours are defined by queen contiguity (cells sharing at least a
#' corner) with row-standardized weights. For offset models the test uses
#' only the cells that entered the fit (non-zero denominator); contiguity is
#' then evaluated among those cells, and cells without any fitted neighbour
#' are dropped from the statistic with a message.
#'
#' A significant result does not invalidate the risk terrain map, but it
#' means reported p-values are anti-conservative and model selection may be
#' too liberal. Recommended follow-up: refit the selected model with a
#' smooth spatial term (e.g. `mgcv::gam(... + s(x, y))`) and report the
#' stability of the RRVs.
#'
#' @param x An `rtm` object.
#' @param type Residual type passed to [stats::residuals.glm()]:
#'   `"pearson"` (default) or `"deviance"`.
#' @param nsim Number of random permutations for the null distribution
#'   (default 999).
#' @param alternative Alternative hypothesis: `"greater"` (positive
#'   autocorrelation, the default and the typical concern), `"less"`, or
#'   `"two.sided"`.
#'
#' @return An object of class `rtm_moran`: a list with `statistic`
#'   (observed Moran's I), `expected` (its null expectation, -1/(n-1)),
#'   `p_value` (permutation p-value), `alternative`, `nsim`, `sims` (the
#'   permuted statistics), `n_cells` (cells used), `n_dropped` (cells
#'   without neighbours), and `type`.
#' @references Moran, P. A. P. (1950). Notes on continuous stochastic
#'   phenomena. *Biometrika*, 37, 17-23.
#'
#'   Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models and
#'   Applications*. London: Pion.
#'
#'   Heffner, J. (2013). Statistics of the RTMDx Utility. In J. Caplan, L.
#'   Kennedy, & E. Piza, *Risk Terrain Modeling Diagnostics Utility User
#'   Manual (Version 1.0)*. Newark, NJ: Rutgers Center on Public Security.
#' @seealso [rtm()], [validate_rtm()]
#' @export
moran_rtm <- function(x, type = c("pearson", "deviance"), nsim = 999,
                      alternative = c("greater", "less", "two.sided")) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  if (nsim < 1) {
    stop("`nsim` must be a positive integer")
  }

  # residuals exist only for the cells that entered the fit
  fit_cells <- if (isTRUE(x$has_offset)) x$grid$offset_count > 0 else rep(TRUE, nrow(x$grid))
  e <- as.numeric(stats::residuals(x$best_model, type = type))
  if (length(e) != sum(fit_cells)) {
    stop(
      "model residuals (", length(e), ") do not match the fitted cells (",
      sum(fit_cells), "); was the object modified?"
    )
  }

  # queen contiguity among the fitted cells (self-intersections removed)
  geom <- sf::st_geometry(x$grid)[fit_cells]
  nb <- sf::st_intersects(geom)
  nb <- lapply(seq_along(nb), function(i) setdiff(nb[[i]], i))

  # drop cells without neighbours; dropping can orphan former neighbours,
  # so repeat until the neighbour graph is stable
  n_dropped <- 0
  repeat {
    island <- lengths(nb) == 0
    if (!any(island)) break
    n_dropped <- n_dropped + sum(island)
    keep <- which(!island)
    reindex <- match(seq_along(nb), keep) # old index -> new index or NA
    nb <- lapply(nb[keep], function(js) {
      js <- reindex[js]
      js[!is.na(js)]
    })
    e <- e[keep]
  }
  if (n_dropped > 0) {
    message(
      n_dropped, " cell(s) without a fitted neighbour dropped from the statistic"
    )
  }

  n <- length(e)
  if (n < 3) {
    stop("fewer than 3 connected cells; Moran's I is not meaningful")
  }

  # sparse edge representation with row-standardized weights: W has S0 = n,
  # so I = z'Wz / z'z with z the centred residuals
  i_idx <- rep(seq_along(nb), lengths(nb))
  j_idx <- unlist(nb)
  w <- 1 / lengths(nb)[i_idx]

  z <- e - mean(e)
  denom <- sum(z^2)
  if (denom == 0) {
    stop("residuals are constant; Moran's I is undefined")
  }

  moran <- function(z) sum(w * z[i_idx] * z[j_idx]) / sum(z^2)

  observed <- moran(z)
  sims <- vapply(seq_len(nsim), function(k) moran(sample(z)), numeric(1))

  p_greater <- (1 + sum(sims >= observed)) / (nsim + 1)
  p_less <- (1 + sum(sims <= observed)) / (nsim + 1)
  p_value <- switch(alternative,
    greater = p_greater,
    less = p_less,
    two.sided = min(1, 2 * min(p_greater, p_less))
  )

  structure(
    list(
      statistic = observed,
      expected = -1 / (n - 1),
      p_value = p_value,
      alternative = alternative,
      nsim = nsim,
      sims = sims,
      n_cells = n,
      n_dropped = n_dropped,
      type = type
    ),
    class = "rtm_moran"
  )
}

#' @export
print.rtm_moran <- function(x, ...) {
  cat("Moran's I test of residual spatial autocorrelation\n")
  cat(
    "Residuals:", x$type, " Cells:", x$n_cells,
    if (x$n_dropped > 0) paste0(" (", x$n_dropped, " without neighbours dropped)"),
    "\n"
  )
  cat(sprintf(
    "Moran's I: %.4f  (expected under independence: %.4f)\n",
    x$statistic, x$expected
  ))
  cat(sprintf(
    "Permutation p-value: %.4g  (alternative: %s, %d permutations)\n",
    x$p_value, x$alternative, x$nsim
  ))
  if (x$alternative == "greater" && x$p_value <= 0.05) {
    cat(
      "\nResiduals are spatially clustered beyond the selected risk factors.\n",
      "Reported p-values are anti-conservative; consider refitting the\n",
      "selected model with a spatial term (e.g. mgcv::gam with s(x, y))\n",
      "and reporting the stability of the RRVs.\n",
      sep = ""
    )
  }
  invisible(x)
}

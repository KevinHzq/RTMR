#' Residual spatial correlogram and cluster-size estimation
#'
#' Computes Moran's I of the model residuals within successive distance
#' bands (a spatial correlogram) to estimate how far residual spatial
#' correlation reaches. The estimated range is the distance at which the
#' band-wise Moran's I decays to the permutation null, and is the
#' data-driven choice for `cluster_size` in [robust_rtm()]: clusters at
#' least as wide as the range contain essentially all of the residual
#' correlation.
#'
#' For each band \eqn{(d_{k-1}, d_k]}, cell pairs whose centroid distance
#' falls in the band define the (row-standardized) spatial weights, and
#' Moran's I is tested against `nsim` random permutations of the residuals
#' (one-sided, positive). The estimated range is the upper edge of the
#' initial consecutive run of significant bands: correlation beyond a gap
#' of non-significant bands is treated as noise rather than an extension of
#' the range. As in [moran_rtm()], offset models use only the cells that
#' entered the fit.
#'
#' By default the bands step by one cell width up to three times the
#' largest candidate spatial influence (at most 15 bands). All cell pairs
#' within the largest band distance are enumerated, so very large grids
#' combined with wide bands can be memory-hungry; reduce `max(breaks)` if
#' needed.
#'
#' @param x An `rtm` object.
#' @param breaks Increasing positive numeric vector of band upper edges, in
#'   CRS units; the bands are `(0, breaks[1]], (breaks[1], breaks[2]], ...`.
#'   Default: one-cell-width steps up to `3 * max(spatial influence)`,
#'   capped at 15 bands.
#' @param nsim Number of residual permutations per band (default 199).
#' @param type Residual type, `"pearson"` (default) or `"deviance"`.
#' @param alpha_level Significance level used to call a band correlated
#'   when estimating the range (default 0.05).
#'
#' @return An object of class `rtm_correlogram`: a list with `bands` (a
#'   data.frame with the band edges, number of cell pairs, Moran's I,
#'   permutation envelope bounds `null_lo`/`null_hi` at `1 - alpha_level`,
#'   and permutation p-value), `range` (estimated correlation range; 0 if
#'   the first band is already non-significant, `Inf` if every band is
#'   significant), `suggested_cluster_size` (equal to `range` when finite
#'   and positive, otherwise `NA`), `expected` (null expectation of I),
#'   `n_cells`, `nsim`, and `alpha_level`.
#' @references Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes:
#'   Models and Applications*. London: Pion.
#' @seealso [moran_rtm()] for the single overall test, [robust_rtm()] for
#'   the cluster-robust inference the result feeds into.
#' @export
correlogram_rtm <- function(x, breaks = NULL, nsim = 199,
                            type = c("pearson", "deviance"),
                            alpha_level = 0.05) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  type <- match.arg(type)
  if (nsim < 1) {
    stop("`nsim` must be a positive integer")
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1 ||
    alpha_level <= 0 || alpha_level >= 1) {
    stop("`alpha_level` must be a single number strictly between 0 and 1")
  }

  fit_cells <- if (isTRUE(x$has_offset)) x$grid$offset_count > 0 else rep(TRUE, nrow(x$grid))
  e <- as.numeric(stats::residuals(x$best_model, type = type))
  if (length(e) != sum(fit_cells)) {
    stop(
      "model residuals (", length(e), ") do not match the fitted cells (",
      sum(fit_cells), "); was the object modified?"
    )
  }

  geom <- sf::st_geometry(x$grid)[fit_cells]
  ctr_sfc <- suppressWarnings(sf::st_centroid(geom))
  ctr <- sf::st_coordinates(ctr_sfc)
  n <- length(e)
  if (n < 3) {
    stop("fewer than 3 fitted cells; a correlogram is not meaningful")
  }

  if (is.null(breaks)) {
    cell_bb <- sf::st_bbox(geom[1])
    cell_width <- as.numeric(cell_bb["xmax"] - cell_bb["xmin"])
    max_d <- 3 * max(x$meta$spatial_influence, na.rm = TRUE)
    if (!is.finite(max_d) || max_d <= 0) {
      max_d <- 10 * cell_width
    }
    step <- max(cell_width, max_d / 15)
    breaks <- seq(step, max_d, by = step)
  }
  breaks <- sort(unique(as.numeric(breaks)))
  if (length(breaks) == 0 || any(!is.finite(breaks)) || breaks[1] <= 0) {
    stop("`breaks` must be positive, finite band upper edges")
  }

  # all directed cell pairs within the widest band, assigned to bands by
  # centroid distance
  hits <- sf::st_is_within_distance(ctr_sfc, dist = max(breaks))
  i_idx <- rep(seq_along(hits), lengths(hits))
  j_idx <- unlist(hits)
  keep <- i_idx != j_idx
  i_idx <- i_idx[keep]
  j_idx <- j_idx[keep]
  d <- sqrt((ctr[i_idx, 1] - ctr[j_idx, 1])^2 + (ctr[i_idx, 2] - ctr[j_idx, 2])^2)
  band <- findInterval(d, c(0, breaks), left.open = TRUE)

  z <- e - mean(e)
  ssz <- sum(z^2)
  if (ssz == 0) {
    stop("residuals are constant; Moran's I is undefined")
  }

  # shared permutations across bands so band-wise nulls are comparable
  perms <- replicate(nsim, sample.int(n))

  lower <- c(0, breaks[-length(breaks)])
  out <- vector("list", length(breaks))
  for (k in seq_along(breaks)) {
    ib <- i_idx[band == k]
    jb <- j_idx[band == k]
    if (length(ib) == 0) {
      out[[k]] <- data.frame(
        lower = lower[k], upper = breaks[k], n_pairs = 0,
        moran = NA_real_, null_lo = NA_real_, null_hi = NA_real_,
        p_value = NA_real_
      )
      next
    }
    deg <- tabulate(ib, nbins = n)
    w <- 1 / deg[ib]
    s0 <- sum(w)
    band_moran <- function(zz) (n / s0) * sum(w * zz[ib] * zz[jb]) / sum(zz^2)

    observed <- band_moran(z)
    sims <- vapply(seq_len(nsim), function(s) band_moran(z[perms[, s]]), numeric(1))
    out[[k]] <- data.frame(
      lower = lower[k], upper = breaks[k], n_pairs = length(ib) / 2,
      moran = observed,
      null_lo = stats::quantile(sims, alpha_level / 2, names = FALSE),
      null_hi = stats::quantile(sims, 1 - alpha_level / 2, names = FALSE),
      p_value = (1 + sum(sims >= observed)) / (nsim + 1)
    )
  }
  bands <- do.call(rbind, out)

  # range: upper edge of the initial consecutive run of significant bands
  sig <- !is.na(bands$p_value) & bands$p_value <= alpha_level & bands$moran > 0
  if (!sig[1]) {
    range_est <- 0
  } else if (all(sig)) {
    range_est <- Inf
  } else {
    range_est <- bands$upper[which(!sig)[1] - 1]
  }

  structure(
    list(
      bands = bands,
      range = range_est,
      suggested_cluster_size = if (is.finite(range_est) && range_est > 0) range_est else NA_real_,
      expected = -1 / (n - 1),
      n_cells = n,
      nsim = nsim,
      alpha_level = alpha_level,
      type = type
    ),
    class = "rtm_correlogram"
  )
}

#' @export
print.rtm_correlogram <- function(x, ...) {
  cat("Residual spatial correlogram (Moran's I by distance band)\n")
  cat(
    "Residuals:", x$type, " Cells:", x$n_cells,
    " Permutations per band:", x$nsim, "\n\n"
  )
  tab <- x$bands
  tab$moran <- round(tab$moran, 4)
  tab$null_lo <- round(tab$null_lo, 4)
  tab$null_hi <- round(tab$null_hi, 4)
  tab$p_value <- signif(tab$p_value, 3)
  print(tab, row.names = FALSE)

  if (x$range == 0) {
    cat(
      "\nNo residual spatial correlation detected in the first band;\n",
      "model-based standard errors are likely adequate.\n",
      sep = ""
    )
  } else if (is.infinite(x$range)) {
    cat(
      "\nResidual correlation persists through the widest band tested;\n",
      "extend `breaks` to find its range, and use a `cluster_size` of at\n",
      "least ", max(x$bands$upper), " in robust_rtm().\n",
      sep = ""
    )
  } else {
    cat(
      "\nEstimated residual correlation range: ", x$range, "\n",
      "Suggested robust_rtm() cluster_size: >= ", x$suggested_cluster_size, "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Plot a residual spatial correlogram
#'
#' @param x An `rtm_correlogram` object from [correlogram_rtm()].
#' @param ... Passed on to [plot()].
#' @export
plot.rtm_correlogram <- function(x, ...) {
  b <- x$bands[!is.na(x$bands$moran), ]
  mid <- (b$lower + b$upper) / 2
  ylim <- range(c(b$moran, b$null_lo, b$null_hi, x$expected), na.rm = TRUE)
  plot(
    mid, b$moran,
    type = "b", pch = 16, ylim = ylim,
    xlab = "Distance band midpoint", ylab = "Moran's I of residuals",
    main = "Residual spatial correlogram", ...
  )
  graphics::polygon(
    c(mid, rev(mid)), c(b$null_lo, rev(b$null_hi)),
    col = grDevices::adjustcolor("grey50", alpha.f = 0.2), border = NA
  )
  graphics::abline(h = x$expected, lty = 2, col = "grey40")
  if (is.finite(x$range) && x$range > 0) {
    graphics::abline(v = x$range, lty = 3, col = "red3")
    graphics::mtext(
      sprintf("estimated range = %g", x$range),
      side = 3, adj = 1, cex = 0.8, col = "red3"
    )
  }
  invisible(x)
}

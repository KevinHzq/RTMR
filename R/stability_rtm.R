#' Bootstrap the selection pipeline: stability and bagged estimates
#'
#' Reruns the full RTMDx-style selection procedure — penalized culling
#' followed by bidirectional stepwise regression under both distributions —
#' on bootstrap resamples of the grid cells, and aggregates which factors
#' are selected and with what coefficients. Unlike sample splitting (see
#' [split_events()]), every resample uses the full data, so nothing is set
#' aside: this is the post-selection tool for sparse settings (few events,
#' rural study areas) where halving the events would starve the selection
#' stage.
#'
#' Three complementary summaries are returned:
#'
#' * **Selection frequency** (Meinshausen & Bühlmann, 2010): the share of
#'   resamples in which each factor is selected at all, and in which each
#'   individual variable (operationalization x spatial influence) is
#'   chosen. A factor selected in 95% of resamples is not a lucky artifact
#'   of one draw; a factor selected in 30% of them should not be presented
#'   as a finding, whatever its p-value.
#' * **Bagged RRV** (Efron, 2014): the factor's coefficient averaged over
#'   *all* resamples, counting 0 when the factor was not selected, then
#'   exponentiated. Averaging over selection shrinks lucky-high estimates
#'   towards 1 in exactly the situations that produce the winner's curse,
#'   giving a selection-aware point estimate without splitting the data.
#' * **Percentile interval**: quantiles of the per-resample coefficient
#'   (including the zeros), reflecting both estimation and selection
#'   uncertainty. With moderate selection frequency the interval rightly
#'   includes 1.
#'
#' Cells can be resampled independently (`resample = "cell"`) or in
#' contiguous spatial blocks (`resample = "block"`, a block bootstrap using
#' square tiles of side `cluster_size`, defaulting to twice the largest
#' candidate spatial influence as in [robust_rtm()]). Block resampling
#' additionally preserves short-range spatial dependence within the
#' resamples and is recommended when [moran_rtm()] flags autocorrelation.
#'
#' A factor's coefficient may be estimated at different spatial influences
#' in different resamples; the aggregation pools them, and the returned
#' variable-level frequencies show how concentrated the choice of
#' influence is. Set a seed before calling for reproducibility. The
#' procedure refits the whole pipeline `nboot` times, so expect a runtime
#' of roughly `nboot` times a single [rtm()] fit (minus the
#' operationalization).
#'
#' @param x An `rtm` object.
#' @param nboot Number of bootstrap resamples (default 200; 50 already
#'   gives a rough stability reading).
#' @param resample `"cell"` (default) resamples grid cells independently;
#'   `"block"` resamples contiguous square tiles of cells.
#' @param cluster_size Tile side length for `resample = "block"`, in CRS
#'   units (default: twice the largest candidate spatial influence).
#' @param nfolds,alpha_level Culling folds and stepwise significance level,
#'   as in [rtm()] (they are not stored in the fitted object; use the same
#'   values as the original call).
#' @param level Level for the percentile intervals (default 0.95).
#' @param verbose Print progress every 10% (default `TRUE`).
#'
#' @return An object of class `rtm_stability`: a list with `factors` (a
#'   data.frame with per-factor selection frequency `freq`, the full-data
#'   `rrv_full` for selected factors, the bagged `rrv_bagged`, and
#'   percentile limits `rrv_lo`/`rrv_hi`), `variables` (per-variable
#'   selection frequencies), `nboot`, `n_failed` (resamples where no model
#'   could be fitted), `resample`, and `level`.
#' @references Efron, B. (2014). Estimation and accuracy after model
#'   selection. *Journal of the American Statistical Association*, 109,
#'   991-1007.
#'
#'   Meinshausen, N., & Bühlmann, P. (2010). Stability selection. *Journal
#'   of the Royal Statistical Society: Series B*, 72, 417-473.
#' @seealso [split_events()] and [refit_rtm()] for the sample-splitting
#'   route when events are plentiful.
#' @export
stability_rtm <- function(x, nboot = 200,
                          resample = c("cell", "block"),
                          cluster_size = NULL,
                          nfolds = 5, alpha_level = 0.05,
                          level = 0.95, verbose = TRUE) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  resample <- match.arg(resample)
  if (nboot < 2) {
    stop("`nboot` must be at least 2")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1")
  }

  fit_cells <- which(
    if (isTRUE(x$has_offset)) x$grid$offset_count > 0 else rep(TRUE, nrow(x$grid))
  )
  data <- x$data[fit_cells, , drop = FALSE]
  meta <- x$meta[, setdiff(names(x$meta), "culled"), drop = FALSE]
  candidates <- meta$variable
  always <- if (!is.null(x$adjustments)) x$adjustments$covariate else character(0)
  offset_col <- if (isTRUE(x$has_offset)) ".log_offset"
  n <- nrow(data)

  # block bootstrap: precompute the tile membership of each fitted cell
  blocks <- NULL
  if (resample == "block") {
    if (is.null(cluster_size)) {
      influences <- x$meta$spatial_influence
      if (all(is.na(influences))) {
        stop("no spatial influence distances found; supply `cluster_size`")
      }
      cluster_size <- 2 * max(influences, na.rm = TRUE)
    }
    ctr <- suppressWarnings(
      sf::st_coordinates(sf::st_centroid(sf::st_geometry(x$grid)))
    )[fit_cells, , drop = FALSE]
    tile <- paste(
      floor((ctr[, 1] - min(ctr[, 1])) / cluster_size),
      floor((ctr[, 2] - min(ctr[, 2])) / cluster_size)
    )
    blocks <- split(seq_len(n), tile)
    if (length(blocks) < 2) {
      stop("block resampling produced fewer than 2 tiles; decrease `cluster_size`")
    }
  }

  coefs <- matrix(
    0,
    nrow = nboot, ncol = length(candidates),
    dimnames = list(NULL, candidates)
  )
  selected_any <- matrix(
    FALSE,
    nrow = nboot, ncol = length(candidates),
    dimnames = list(NULL, candidates)
  )
  failed <- logical(nboot)

  for (b in seq_len(nboot)) {
    idx <- if (resample == "cell") {
      sample.int(n, n, replace = TRUE)
    } else {
      unlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)],
        use.names = FALSE
      )
    }
    df <- data[idx, , drop = FALSE]

    # variables constant in this resample cannot enter
    usable <- candidates[
      vapply(candidates, function(v) length(unique(df[[v]])) > 1, logical(1))
    ]

    kept <- usable
    if (length(usable) > 1) {
      kept <- tryCatch(
        cull_variables(
          df[usable], df$outcome_count,
          nfolds = nfolds,
          offset = if (!is.null(offset_col)) df[[offset_col]],
          adjust = if (length(always) > 0) df[always]
        ),
        error = function(e) usable
      )
    }

    results <- list()
    for (fam in c("poisson", "nb")) {
      results[[fam]] <- rtm_stepwise(
        data = df, outcome = "outcome_count",
        meta = meta[meta$variable %in% kept, , drop = FALSE],
        family = fam, alpha_level = alpha_level,
        offset_col = offset_col, always = always
      )
    }
    results <- Filter(Negate(is.null), results)
    if (length(results) == 0) {
      failed[b] <- TRUE
      next
    }
    best <- results[[which.min(vapply(results, `[[`, numeric(1), "bic"))]]
    if (length(best$variables) > 0) {
      est <- stats::coef(best$fit)[best$variables]
      coefs[b, best$variables] <- est
      selected_any[b, best$variables] <- TRUE
    }
    if (verbose && b %% max(1, nboot %/% 10) == 0) {
      message("  bootstrap ", b, "/", nboot)
    }
  }

  ok <- !failed
  n_ok <- sum(ok)
  if (n_ok < 2) {
    stop("fewer than 2 bootstrap resamples produced a model; check the fit")
  }

  variables <- data.frame(
    meta[, c("variable", "factor", "operation", "spatial_influence")],
    freq = colMeans(selected_any[ok, , drop = FALSE]),
    row.names = NULL
  )

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  fac_names <- unique(meta$factor)
  factors <- do.call(rbind, lapply(fac_names, function(f) {
    vs <- meta$variable[meta$factor == f]
    fac_coef <- rowSums(coefs[ok, vs, drop = FALSE])
    fac_sel <- rowSums(selected_any[ok, vs, drop = FALSE]) > 0
    q <- stats::quantile(fac_coef, probs, names = FALSE)
    data.frame(
      factor = f,
      freq = mean(fac_sel),
      rrv_full = if (f %in% x$selected$factor) {
        x$selected$rrv[match(f, x$selected$factor)]
      } else {
        NA_real_
      },
      rrv_bagged = exp(mean(fac_coef)),
      rrv_lo = exp(q[1]),
      rrv_hi = exp(q[2])
    )
  }))
  factors <- factors[order(-factors$freq), ]
  rownames(factors) <- NULL

  structure(
    list(
      factors = factors,
      variables = variables,
      nboot = nboot,
      n_failed = sum(failed),
      resample = resample,
      cluster_size = cluster_size,
      level = level
    ),
    class = "rtm_stability"
  )
}

#' @export
print.rtm_stability <- function(x, ...) {
  cat("Bootstrap stability of the risk terrain selection pipeline\n")
  cat(
    "Resamples:", x$nboot - x$n_failed, "successful of", x$nboot,
    sprintf(
      " (%s bootstrap%s)\n\n", x$resample,
      if (x$resample == "block") sprintf(", tile %g", x$cluster_size) else ""
    )
  )

  tab <- x$factors
  tab$freq <- round(tab$freq, 2)
  tab$rrv_full <- round(tab$rrv_full, 3)
  tab$rrv_bagged <- round(tab$rrv_bagged, 3)
  tab$rrv_lo <- round(tab$rrv_lo, 3)
  tab$rrv_hi <- round(tab$rrv_hi, 3)
  print(tab, row.names = FALSE)

  cat(
    "\nfreq: share of resamples selecting the factor. rrv_bagged averages\n",
    "the coefficient over all resamples (0 when unselected), shrinking\n",
    "selection-inflated estimates; the interval includes selection\n",
    "uncertainty. Factors with low freq should not be reported as findings.\n",
    sep = ""
  )
  invisible(x)
}

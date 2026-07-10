#' Randomly split events for selection-clean inference
#'
#' Splits an event point layer into a training and a test layer by
#' independently assigning each event to the training set with probability
#' `prop` (Poisson thinning, also known as count splitting; Neufeld et al.,
#' 2024). Select the model on the training half with [rtm()], then
#' re-estimate the selected specification on the test half with
#' [refit_rtm()]: because thinning a Poisson process yields two independent
#' processes over the *same* cells and time period, the test-half
#' coefficients are untainted by model selection and unaffected by temporal
#' change in the effects.
#'
#' Under Poisson variation the two halves are exactly independent. Under
#' overdispersion (the negative binomial case, e.g. events that trigger
#' nearby events) they are mildly positively dependent, making test-half
#' inference slightly optimistic — still far closer to honest than reusing
#' the selection data.
#'
#' Splitting halves the events available to each stage, so expect less
#' power in both selection and inference; `prop` can be raised above 0.5 to
#' favour selection or lowered to favour inference.
#'
#' @param pt An sf or sfc object with the event points.
#' @param prop Probability that an event is assigned to the training half
#'   (default 0.5).
#'
#' @return A list with sf/sfc elements `train` and `test`, a partition of
#'   the rows of `pt`.
#' @references Neufeld, A., Dharamshi, A., Gao, L. L., & Witten, D. (2024).
#'   Data thinning to avoid double dipping. *Journal of the American
#'   Statistical Association*, 119, 2416-2427.
#' @seealso [refit_rtm()]
#' @export
split_events <- function(pt, prop = 0.5) {
  if (!inherits(pt, c("sf", "sfc"))) {
    stop("`pt` must be an sf or sfc object")
  }
  if (!is.numeric(prop) || length(prop) != 1 || prop <= 0 || prop >= 1) {
    stop("`prop` must be a single number strictly between 0 and 1")
  }
  n <- length(sf::st_geometry(pt))
  if (n < 2) {
    stop("`pt` must contain at least 2 events to split")
  }
  in_train <- stats::runif(n) < prop
  if (inherits(pt, "sf")) {
    list(train = pt[in_train, ], test = pt[!in_train, ])
  } else {
    list(train = pt[in_train], test = pt[!in_train])
  }
}

#' Re-estimate a fitted risk terrain model on new events
#'
#' Takes the *frozen* specification of a fitted [rtm()] model — the
#' selected variables, adjustment covariates, offset structure, and
#' distribution family — and re-estimates its coefficients on a new set of
#' outcome events, without repeating any selection. Because the new data
#' played no part in choosing the specification, the refitted p-values and
#' confidence intervals are ordinary ones, free of the post-selection
#' caveat attached to the original fit.
#'
#' What a difference between original and refitted RRVs *means* depends on
#' where the new events come from:
#'
#' * **Same period, thinned split** (see [split_events()]): temporal change
#'   is excluded by construction, so systematic attenuation of the
#'   refitted RRVs towards 1 is the signature of selection bias (winner's
#'   curse) in the original estimates.
#' * **A later period** (temporal holdout): differences combine selection
#'   bias *and* genuine change in the effects over time. This is the
#'   relevant quantity for prospective use ("do these effects hold now?"),
#'   but it cannot separate the two sources; running both designs can —
#'   attenuation already present in the thinned split is selection, any
#'   additional change in the temporal refit is drift.
#'
#' The refit keeps the original grid, operationalized risk factor
#' variables, and cell-level covariates: only the outcome counts (and
#' optionally the offset denominator) are replaced. Risk factor layers are
#' therefore assumed unchanged; if the environment itself changed between
#' periods (new bars, closed clinics), rebuild the model instead. The
#' distribution family is kept fixed; for a negative binomial model the
#' dispersion parameter is re-estimated on the new data.
#'
#' @param x An `rtm` object.
#' @param new_events An sf or sfc object with the new outcome event points,
#'   in the same CRS as the model grid.
#' @param new_offset Replacement denominator for offset models: an sf/sfc
#'   point layer (counted per cell) or a numeric vector with one value per
#'   grid cell. If `NULL` (default) the original denominator is reused,
#'   which is appropriate for same-period thinned splits where only the
#'   outcome events were split.
#' @param level Confidence level for the refitted RRV intervals (default
#'   0.95).
#'
#' @return An object of class `rtm_refit`: a list with `selected` (the
#'   selected-variables table comparing `rrv_original` with `rrv_refit`
#'   and its `ci_lower`/`ci_upper` and `p_refit`), `adjustments` (original
#'   vs refitted coefficients of adjustment covariates, `NULL` if none),
#'   `intercept` (refitted), `family`, `n_events`, `n_cells` (cells used
#'   in the refit), `converged`, and `level`.
#' @seealso [split_events()] to create a same-period split,
#'   [validate_rtm()] for the discrimination side of holdout checking.
#' @export
refit_rtm <- function(x, new_events, new_offset = NULL, level = 0.95) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1")
  }
  if (nrow(x$selected) == 0) {
    stop("the model has no selected risk factors to re-estimate")
  }

  counts <- count_points(new_events, x$grid)
  if (sum(counts) == 0) {
    stop("no new events fall within the model grid; check the CRS and study area")
  }

  data <- x$data
  data$outcome_count <- counts

  offset_col <- NULL
  if (isTRUE(x$has_offset)) {
    offset_col <- ".log_offset"
    if (is.null(new_offset)) {
      offset_count <- x$grid$offset_count
    } else {
      offset_count <- if (inherits(new_offset, c("sf", "sfc"))) {
        count_points(new_offset, x$grid)
      } else {
        as.numeric(new_offset)
      }
      if (length(offset_count) != nrow(x$grid)) {
        stop(
          "`new_offset` must have one value per grid cell (",
          nrow(x$grid), " cells)"
        )
      }
      if (anyNA(offset_count) || any(offset_count < 0) || all(offset_count == 0)) {
        stop("`new_offset` counts must be non-negative, non-missing, and not all zero")
      }
      data$.log_offset <- ifelse(offset_count > 0, log(offset_count), 0)
    }
    fit_cells <- offset_count > 0
    dropped <- sum(counts[!fit_cells])
    if (dropped > 0) {
      warning(
        dropped, " new event(s) fall in cells with a zero denominator ",
        "and are excluded from the refit"
      )
    }
  } else {
    if (!is.null(new_offset)) {
      stop("`new_offset` only applies to models fitted with an offset")
    }
    fit_cells <- rep(TRUE, nrow(x$grid))
  }

  always <- if (!is.null(x$adjustments)) x$adjustments$covariate else character(0)
  variables <- x$selected$variable

  fit <- fit_count_model(
    c(always, variables), data[fit_cells, ],
    outcome = "outcome_count", family = x$family, offset_col = offset_col
  )
  if (is.null(fit)) {
    stop("the selected specification could not be refitted on the new events")
  }

  est <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  zq <- stats::qnorm(1 - (1 - level) / 2)

  term_table <- function(terms) {
    idx <- match(terms, names(est))
    data.frame(
      coef_refit = est[idx],
      se_refit = se[idx],
      p_refit = 2 * stats::pnorm(-abs(est[idx] / se[idx])),
      ci_lower = est[idx] - zq * se[idx],
      ci_upper = est[idx] + zq * se[idx],
      row.names = NULL
    )
  }

  tab <- term_table(variables)
  selected <- x$selected[, c("factor", "type", "operation", "spatial_influence", "variable")]
  selected$rrv_original <- x$selected$rrv
  selected$rrv_refit <- exp(tab$coef_refit)
  selected$ci_lower <- exp(tab$ci_lower)
  selected$ci_upper <- exp(tab$ci_upper)
  selected$p_refit <- tab$p_refit

  adjustments <- NULL
  if (length(always) > 0) {
    adj <- term_table(always)
    adjustments <- data.frame(
      covariate = always,
      coef_original = x$adjustments$coefficient,
      coef_refit = adj$coef_refit,
      se_refit = adj$se_refit,
      p_refit = adj$p_refit,
      row.names = NULL
    )
  }

  structure(
    list(
      selected = selected,
      adjustments = adjustments,
      intercept = unname(est["(Intercept)"]),
      family = x$family,
      n_events = sum(counts[fit_cells]),
      n_cells = sum(fit_cells),
      converged = isTRUE(fit$converged),
      level = level
    ),
    class = "rtm_refit"
  )
}

#' @export
print.rtm_refit <- function(x, ...) {
  cat("Re-estimated risk terrain model (frozen specification)\n")
  cat(
    "Family:", if (x$family == "nb") "negative binomial" else "Poisson",
    " New events:", x$n_events, " Cells:", x$n_cells, "\n"
  )
  if (!x$converged) {
    cat("Warning: the refit did not fully converge; interpret with caution\n")
  }
  cat("\nSelected risk factors, original vs refitted:\n")
  tab <- x$selected[, c(
    "factor", "operation", "spatial_influence",
    "rrv_original", "rrv_refit", "ci_lower", "ci_upper", "p_refit"
  )]
  tab$rrv_original <- round(tab$rrv_original, 3)
  tab$rrv_refit <- round(tab$rrv_refit, 3)
  tab$ci_lower <- round(tab$ci_lower, 3)
  tab$ci_upper <- round(tab$ci_upper, 3)
  tab$p_refit <- signif(tab$p_refit, 3)
  print(tab, row.names = FALSE)

  if (!is.null(x$adjustments)) {
    cat("\nAdjustment covariates:\n")
    adj <- x$adjustments
    adj$coef_original <- round(adj$coef_original, 4)
    adj$coef_refit <- round(adj$coef_refit, 4)
    adj$se_refit <- round(adj$se_refit, 4)
    adj$p_refit <- signif(adj$p_refit, 3)
    print(adj, row.names = FALSE)
  }

  cat(
    "\nRefitted p-values and intervals are free of selection effects.\n",
    "Interpreting original-vs-refit differences depends on the data source\n",
    "(same-period split: selection bias; later period: selection + change\n",
    "over time) - see ?refit_rtm.\n",
    sep = ""
  )
  invisible(x)
}

#' Validate a risk terrain model against holdout events
#'
#' Assesses how well a fitted [rtm()] model predicts events it has not seen,
#' typically from a later time period (fit on period 1, validate on period
#' 2). Two complementary views of predictive accuracy are reported:
#'
#' * **Ranked capture / PAI.** Cells are ranked by relative risk score and
#'   the share of holdout events falling in the riskiest `area_pct` percent
#'   of cells is computed. The Predictive Accuracy Index (PAI; Chainey,
#'   Tompson & Uhlig, 2008) is the percent of events captured divided by the
#'   percent of area flagged, so PAI = 8 means the flagged cells catch events
#'   at 8 times the study-area-wide rate.
#' * **ROC / AUC.** Each cell is treated as a binary case (at least one
#'   holdout event vs. none) and the risk score as the classifier. The AUC
#'   summarizes discrimination: 0.5 is chance, 1 is perfect.
#'
#' The observed/expected ratio compares the total holdout events with the
#' model's total predicted count; it is only meaningful as a calibration
#' check when the holdout period has the same duration as the fitting period
#' (the model's predictions represent the expectation for an equal-length
#' period).
#'
#' For offset (rate) models, cells are ranked by their relative rate
#' multiplier, and plain capture/PAI compare against *area*, not against the
#' denominator: a high-lethality cell with few denominator events may
#' rightly capture few holdout events. Supplying `new_offset` (the holdout
#' period's denominator, e.g. all overdose events in period 2) makes the
#' validation denominator-aware: the capture table gains the share of the
#' holdout *denominator* in the flagged cells and an adjusted PAI
#' (`adj_pai`), the ratio of the percent of outcome events captured to the
#' percent of denominator events in the same cells. `adj_pai` above 1 means
#' outcome events concentrate in the cells ranked most risky *beyond* what
#' the denominator distribution alone would produce — for a case-fatality
#' model, that events occurring in flagged cells are disproportionately
#' fatal. The observed/expected ratio is then also computed against the
#' model's expected counts given the holdout denominator (a genuine
#' calibration check that no longer depends on equal period durations).
#'
#' @param x An `rtm` object.
#' @param new_events An sf or sfc object with holdout event points, in the
#'   same CRS as the model grid.
#' @param area_pct Percentages of the study area (riskiest cells first) at
#'   which to evaluate event capture.
#' @param new_offset Optional holdout-period denominator for models fitted
#'   with an offset: an sf/sfc point layer (counted per cell) or a numeric
#'   vector with one value per grid cell.
#'
#' @return An object of class `rtm_validation`: a list with elements
#'   `capture` (data.frame with `area_pct`, `n_cells`, `n_events`,
#'   `pct_events`, and `pai`; with `new_offset` also `n_denom`, `pct_denom`,
#'   and `adj_pai`), `auc`, `roc` (data.frame with `fpr` and `tpr`
#'   for plotting), `n_events`, `n_cells`, `observed_expected` (ratio of
#'   total holdout events to total expected events), and `counts` (holdout
#'   events per cell, in grid order).
#' @references Chainey, S., Tompson, L., & Uhlig, S. (2008). The utility of
#'   hotspot mapping for predicting spatial patterns of crime. *Security
#'   Journal*, 21, 4-28.
#' @export
validate_rtm <- function(x, new_events, area_pct = c(1, 5, 10, 20),
                         new_offset = NULL) {
  score <- x$grid$relrisk
  counts <- count_points(new_events, x$grid)

  n_cells <- length(score)
  n_events <- sum(counts)
  if (n_events == 0) {
    stop("no holdout events fall within the model grid; check the CRS and study area")
  }

  # holdout denominator and model-expected counts given that denominator
  denom <- NULL
  expected <- x$grid$prediction
  if (!is.null(new_offset)) {
    if (!isTRUE(x$has_offset)) {
      stop("`new_offset` only applies to models fitted with an offset")
    }
    denom <- if (inherits(new_offset, c("sf", "sfc"))) {
      count_points(new_offset, x$grid)
    } else {
      as.numeric(new_offset)
    }
    if (length(denom) != n_cells) {
      stop("`new_offset` must have one value per grid cell (", n_cells, " cells)")
    }
    if (sum(denom) == 0) {
      stop("no holdout denominator events fall within the model grid")
    }
    newdata <- x$data
    newdata$.log_offset <- ifelse(denom > 0, log(denom), 0)
    expected <- ifelse(
      denom > 0,
      as.numeric(stats::predict(x$best_model, newdata = newdata, type = "response")),
      0
    )
  }

  # ranked capture and PAI
  ord <- order(-score)
  cum_events <- cumsum(counts[ord])
  cum_denom <- if (!is.null(denom)) cumsum(denom[ord])
  capture <- do.call(rbind, lapply(area_pct, function(p) {
    n_top <- max(1L, ceiling(p / 100 * n_cells))
    pct_area <- 100 * n_top / n_cells
    pct_events <- 100 * cum_events[n_top] / n_events
    out <- data.frame(
      area_pct = p,
      n_cells = n_top,
      n_events = cum_events[n_top],
      pct_events = pct_events,
      pai = pct_events / pct_area
    )
    if (!is.null(denom)) {
      pct_denom <- 100 * cum_denom[n_top] / sum(denom)
      out$n_denom <- cum_denom[n_top]
      out$pct_denom <- pct_denom
      out$adj_pai <- ifelse(pct_denom > 0, pct_events / pct_denom, NA_real_)
    }
    out
  }))

  # ROC over cells: positive = cell received at least one holdout event
  positive <- counts > 0
  n_pos <- sum(positive)
  n_neg <- n_cells - n_pos

  if (n_pos == 0 || n_neg == 0) {
    warning("all cells are of one class; AUC is undefined")
    auc <- NA_real_
    roc <- data.frame(fpr = c(0, 1), tpr = c(0, 1))
  } else {
    # Mann-Whitney formulation; midranks handle tied scores
    r <- rank(score)
    auc <- (sum(r[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    # ROC curve: sweep the threshold down through distinct score values
    tp <- cumsum(positive[ord])
    fp <- cumsum(!positive[ord])
    last_of_tie <- c(diff(score[ord]) != 0, TRUE)
    roc <- data.frame(
      fpr = c(0, fp[last_of_tie] / n_neg),
      tpr = c(0, tp[last_of_tie] / n_pos)
    )
  }

  structure(
    list(
      capture = capture,
      auc = auc,
      roc = roc,
      n_events = n_events,
      n_cells = n_cells,
      observed_expected = n_events / sum(expected),
      denominator_aware = !is.null(denom),
      counts = counts
    ),
    class = "rtm_validation"
  )
}

#' @export
print.rtm_validation <- function(x, ...) {
  cat("Risk terrain model validation\n")
  cat(
    "Holdout events:", x$n_events, " Cells:", x$n_cells,
    sprintf(" Observed/expected: %.2f\n\n", x$observed_expected)
  )

  cat("Event capture by riskiest cells:\n")
  tab <- x$capture
  tab$pct_events <- round(tab$pct_events, 1)
  tab$pai <- round(tab$pai, 2)
  if (isTRUE(x$denominator_aware)) {
    tab$pct_denom <- round(tab$pct_denom, 1)
    tab$adj_pai <- round(tab$adj_pai, 2)
  }
  print(tab, row.names = FALSE)
  if (isTRUE(x$denominator_aware)) {
    cat("adj_pai: outcome capture relative to denominator share (>1 = flagged cells are disproportionately risky)\n")
  }

  cat(sprintf("\nAUC: %.3f\n", x$auc))
  invisible(x)
}

#' Plot the ROC curve of a risk terrain model validation
#'
#' @param x An `rtm_validation` object from [validate_rtm()].
#' @param ... Passed on to [plot()].
#' @export
plot.rtm_validation <- function(x, ...) {
  plot(
    x$roc$fpr, x$roc$tpr,
    type = "l",
    xlab = "False positive rate (share of no-event cells flagged)",
    ylab = "True positive rate (share of event cells flagged)",
    main = sprintf("ROC curve (AUC = %.3f)", x$auc),
    ...
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")
  invisible(x)
}

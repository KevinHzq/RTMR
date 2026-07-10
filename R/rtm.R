#' Fit a risk terrain model
#'
#' Runs the full risk terrain modeling pipeline, replicating the statistical
#' procedure of the RTMDx Utility (Heffner, 2013):
#'
#' 1. A grid of cells is created over the study area and the outcome events
#'    are counted per cell.
#' 2. Each risk factor is operationalized into candidate binary variables at
#'    a range of spatial influences (proximity thresholds and/or high-density
#'    areas from Epanechnikov kernel densities), see [operationalize()].
#' 3. The candidate set is culled with a cross-validated, sign-constrained
#'    penalized Poisson regression to guard against spurious correlations,
#'    see [cull_variables()].
#' 4. Bidirectional stepwise regression by BIC selects the final risk factors
#'    and their optimal spatial influences, run under both a Poisson and a
#'    negative binomial model (to allow for overdispersion); the distribution
#'    with the best BIC wins. See [rtm_stepwise()].
#' 5. The final model's fitted values are the expected event counts per cell;
#'    dividing by the smallest expectation gives the relative risk score
#'    (starting at 1 = least risky cell). Each selected factor's exponentiated
#'    coefficient is its relative risk value (RRV).
#'
#' @param outcome An sf or sfc object with the outcome event points (e.g.
#'   crimes, or fatal overdoses).
#' @param offset Optional denominator turning the model into a rate model:
#'   either an sf/sfc point layer whose per-cell counts form the denominator
#'   (e.g. all overdose events when the outcome is fatal overdoses, giving a
#'   case-fatality model), or a numeric vector with one value per grid cell
#'   (e.g. resident population). Internally entered as a log offset in the
#'   regressions, so RRVs become adjusted rate ratios (e.g. case-fatality
#'   ratios) and `relrisk` becomes the relative rate multiplier surface.
#'   Cells with a zero denominator are excluded from model fitting (their
#'   rate is undefined), but `relrisk` is still computed for them from their
#'   risk factor exposures.
#' @param factors A named list of risk factors. Each element is either an sf
#'   or sfc object, or a list with element `data` (the sf object) plus
#'   optional per-factor overrides of `operation`, `max_blocks`, `increment`,
#'   and `type` (e.g. `type = "protective"`). Point layers support any
#'   operation; line and polygon layers support `operation = "proximity"`
#'   only (see [operationalize()]).
#' @param covariates Optional named list of pre-computed cell-level
#'   covariates (e.g. raster values aggregated to the grid). Each element
#'   must be a list with `values` (a numeric vector, one value per grid
#'   cell) and `role`:
#'   \describe{
#'     \item{`role = "candidate"`}{The covariate competes for selection like
#'       an operationalized risk factor and, if selected, is reported with
#'       an RRV. Give it a `type` (`"aggravating"`, the default, or
#'       `"protective"`; protective values are negated internally as in
#'       [operationalize()]). Binary 0/1 values are recommended so the RRV
#'       reads as an exposed/unexposed rate ratio; continuous values are
#'       allowed with a warning, and their RRV is a rate ratio per one-unit
#'       increase.}
#'     \item{`role = "adjustment"`}{The covariate is included in every
#'       candidate model (a confounder held constant during selection, e.g.
#'       area-level deprivation). It is unpenalized in the culling step,
#'       may take either coefficient sign, is not subject to the
#'       significance check, and is reported separately from the RRV table.
#'       Continuous values are fine here. Categorical values (factor or
#'       character, e.g. land-use class per cell) are also supported: they
#'       are dummy-coded against the level named in `reference` (default:
#'       the first level), and each non-reference level's coefficient is
#'       the adjusted log rate ratio versus the reference. Categorical
#'       values are not allowed for candidates, which need a directional
#'       exposure; hypothesize a specific class with a binary indicator
#'       instead, e.g. `values = as.integer(landuse == "commercial")`.}
#'   }
#'   Covariate values must be complete (no `NA`).
#' @param boundary An sf object with the study area polygon(s). All layers
#'   must share a projected CRS.
#' @param cell_size Grid cell size in CRS units. RTMDx recommends half the
#'   block length.
#' @param block_length Average block length in CRS units.
#' @param operation,max_blocks,increment,type Defaults applied to every
#'   factor unless overridden per factor; see [operationalize()].
#' @param cull If `TRUE` (default), run the penalized-regression culling step.
#' @param nfolds Cross-validation folds for culling (default 5).
#' @param alpha_level Significance level for the stepwise validity check.
#' @param verbose Print progress messages (default `TRUE`).
#' @param ... Passed on to [create_grid()] (and hence [sf::st_make_grid()]),
#'   e.g. `square = FALSE`.
#'
#' @return An object of class `rtm`, a list with:
#'   \describe{
#'     \item{best_model}{The final fitted model (`glm` or `negbin`).}
#'     \item{family}{`"poisson"` or `"nb"`.}
#'     \item{bic}{BIC of the final model.}
#'     \item{selected}{data.frame of selected variables with their factor,
#'       operation, spatial influence, coefficient, p-value, and relative
#'       risk value. Selected candidate covariates appear with operation
#'       `"covariate"`.}
#'     \item{adjustments}{data.frame with the coefficient and p-value of
#'       each adjustment covariate (`NULL` if none).}
#'     \item{grid}{The sf grid with columns `outcome_count`, `prediction`
#'       (expected event count), and `relrisk` (relative risk score). For
#'       offset models also `offset_count`, with `prediction` the expected
#'       outcome count given the denominator and `relrisk` the relative rate
#'       (e.g. relative case fatality), rescaled to start at 1.}
#'     \item{intercept}{Model intercept.}
#'     \item{data}{The cell-by-variable model data.}
#'     \item{meta}{Metadata for all candidate variables, with a `culled`
#'       column.}
#'   }
#' @references Heffner, J. (2013). Statistics of the RTMDx Utility. In J.
#'   Caplan, L. Kennedy, & E. Piza, *Risk Terrain Modeling Diagnostics
#'   Utility User Manual (Version 1.0)*. Newark, NJ: Rutgers Center on Public
#'   Security.
#' @export
rtm <- function(outcome, factors, boundary,
                cell_size, block_length,
                offset = NULL,
                covariates = NULL,
                operation = c("proximity", "density", "both"),
                max_blocks = 3,
                increment = c("whole", "half"),
                type = c("aggravating", "protective"),
                cull = TRUE,
                nfolds = 5,
                alpha_level = 0.05,
                verbose = TRUE,
                ...) {
  operation <- match.arg(operation)
  increment <- match.arg(increment)
  type <- match.arg(type)

  if (is.null(names(factors)) || any(!nzchar(names(factors)))) {
    stop("`factors` must be a named list of risk factor layers")
  }

  say <- function(...) if (verbose) message(...)

  # 1. grid and outcome counts
  say("Creating grid and counting outcome events...")
  grid <- create_grid(boundary, cellsize = cell_size, ...)
  outcome_count <- count_points(outcome, grid)

  # optional denominator: cells with zero denominator have an undefined rate
  # and are excluded from fitting
  if (!is.null(offset)) {
    offset_count <- if (inherits(offset, c("sf", "sfc"))) {
      count_points(offset, grid)
    } else {
      as.numeric(offset)
    }
    if (length(offset_count) != nrow(grid)) {
      stop("`offset` must have one value per grid cell (", nrow(grid), " cells)")
    }
    if (any(offset_count < 0) || all(offset_count == 0)) {
      stop("`offset` counts must be non-negative and not all zero")
    }
    fit_cells <- offset_count > 0
    say(
      sum(fit_cells), " of ", nrow(grid),
      " cells have a non-zero denominator and enter the model"
    )
    dropped_events <- sum(outcome_count[!fit_cells])
    if (dropped_events > 0) {
      warning(
        dropped_events, " outcome event(s) fall in cells with a zero ",
        "denominator and are excluded from fitting; check that the ",
        "denominator layer covers the outcome events"
      )
    }
  } else {
    offset_count <- NULL
    fit_cells <- rep(TRUE, nrow(grid))
  }

  # 2. operationalize risk factors into candidate binary variables
  say("Operationalizing risk factors...")
  op_list <- lapply(names(factors), function(nm) {
    f <- factors[[nm]]
    if (!inherits(f, c("sf", "sfc"))) {
      spec <- f
      f <- spec$data
    } else {
      spec <- list()
    }
    operationalize(
      pt = f, grid = grid, name = nm,
      block_length = block_length,
      max_blocks = spec$max_blocks %||% max_blocks,
      increment = spec$increment %||% increment,
      operation = spec$operation %||% operation,
      type = spec$type %||% type
    )
  })

  meta <- do.call(rbind, lapply(op_list, attr, "rtm_meta"))
  x <- do.call(cbind, op_list)

  # pre-computed cell-level covariates: candidates join the selection pool,
  # adjustments are carried separately and included in every model
  adjust <- NULL
  if (!is.null(covariates)) {
    if (is.null(names(covariates)) || any(!nzchar(names(covariates)))) {
      stop("`covariates` must be a named list")
    }
    for (nm in names(covariates)) {
      spec <- covariates[[nm]]
      if (!is.list(spec) || is.null(spec$values) || is.null(spec$role)) {
        stop(
          "covariate `", nm, "` must be a list with `values` and ",
          "`role` (\"candidate\" or \"adjustment\")"
        )
      }
      role <- match.arg(spec$role, c("candidate", "adjustment"))
      vals <- spec$values
      if (length(vals) != nrow(grid)) {
        stop(
          "covariate `", nm, "` must have one value per grid cell (",
          nrow(grid), " cells)"
        )
      }
      if (anyNA(vals)) {
        stop(
          "covariate `", nm, "` contains missing values; ",
          "fill or drop them before modeling"
        )
      }
      categorical <- is.character(vals) || is.factor(vals)
      if (!categorical && !is.null(spec$reference)) {
        warning("`reference` for covariate `", nm, "` is ignored: it only applies to categorical values")
      }
      v <- make.names(nm)
      if (role == "candidate") {
        if (categorical) {
          stop(
            "candidate covariate `", nm, "` is categorical; candidates need a ",
            "directional exposure - pass a binary indicator for the class of ",
            "interest instead, e.g. values = as.integer(", nm, " == \"<level>\")"
          )
        }
        vals <- as.numeric(vals)
        ctype <- match.arg(spec$type %||% "aggravating", c("aggravating", "protective"))
        if (!all(vals %in% c(0, 1))) {
          warning(
            "candidate covariate `", nm, "` is not binary 0/1; ",
            "its RRV will be a rate ratio per one-unit increase"
          )
        }
        x[[v]] <- if (ctype == "protective") -vals else vals
        meta <- rbind(meta, data.frame(
          variable = v, factor = v, operation = "covariate",
          spatial_influence = NA_real_, type = ctype,
          stringsAsFactors = FALSE
        ))
      } else if (categorical) {
        f <- droplevels(factor(vals))
        if (nlevels(f) < 2) {
          warning("dropping constant adjustment covariate: ", nm)
          next
        }
        ref <- spec$reference %||% levels(f)[1]
        if (!ref %in% levels(f)) {
          stop(
            "reference level \"", ref, "\" not found in covariate `", nm,
            "` (levels: ", paste(levels(f), collapse = ", "), ")"
          )
        }
        f <- stats::relevel(f, ref = ref)
        dummies <- stats::model.matrix(~f)[, -1, drop = FALSE]
        colnames(dummies) <- make.names(paste0(nm, "_", sub("^f", "", colnames(dummies))))
        for (dv in colnames(dummies)) {
          adjust[[dv]] <- dummies[, dv]
        }
      } else {
        adjust[[v]] <- as.numeric(vals)
      }
    }
    if (!is.null(adjust)) {
      adjust <- as.data.frame(adjust)
      adj_constant <- vapply(adjust, function(v) length(unique(v[fit_cells])) < 2, logical(1))
      if (any(adj_constant)) {
        warning(
          "dropping constant adjustment covariate(s): ",
          paste(names(adjust)[adj_constant], collapse = ", ")
        )
        adjust <- adjust[!adj_constant]
        if (ncol(adjust) == 0) adjust <- NULL
      }
    }
  }

  # drop variables constant across the fitted cells (e.g. a spatial
  # influence covering every cell)
  constant <- vapply(x, function(v) length(unique(v[fit_cells])) < 2, logical(1))
  if (any(constant)) {
    say(
      "Dropping ", sum(constant), " constant variable(s): ",
      paste(names(x)[constant], collapse = ", ")
    )
    x <- x[!constant]
    meta <- meta[!constant, ]
  }
  say(nrow(meta), " candidate variables for ", length(factors), " risk factors")

  # 3. cull variables with cross-validated penalized Poisson regression
  if (cull && ncol(x) > 1) {
    say("Culling variables with penalized regression...")
    kept <- cull_variables(
      x[fit_cells, ], outcome_count[fit_cells],
      nfolds = nfolds,
      offset = if (!is.null(offset_count)) log(offset_count[fit_cells]),
      adjust = if (!is.null(adjust)) adjust[fit_cells, , drop = FALSE]
    )
    say(length(kept), " variable(s) survived culling")
  } else {
    kept <- names(x)
  }
  meta$culled <- !(meta$variable %in% kept)

  model_data <- cbind(outcome_count = outcome_count, x)
  always <- character(0)
  if (!is.null(adjust)) {
    model_data <- cbind(model_data, adjust)
    always <- names(adjust)
  }
  offset_col <- NULL
  if (!is.null(offset_count)) {
    # log offset for the fitted cells; 0 elsewhere so predictions on the
    # full grid yield the rate multiplier
    log_off <- numeric(nrow(grid))
    log_off[fit_cells] <- log(offset_count[fit_cells])
    model_data$.log_offset <- log_off
    offset_col <- ".log_offset"
  }

  # 4. bidirectional stepwise regression by BIC, Poisson vs negative binomial
  results <- list()
  for (fam in c("poisson", "nb")) {
    say("Stepwise regression (", fam, ")...")
    results[[fam]] <- rtm_stepwise(
      data = model_data[fit_cells, ], outcome = "outcome_count",
      meta = meta[!meta$culled, ],
      family = fam, alpha_level = alpha_level,
      offset_col = offset_col, always = always, verbose = verbose
    )
  }
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    stop("no model could be fitted; check the outcome counts")
  }

  best <- results[[which.min(vapply(results, `[[`, numeric(1), "bic"))]]
  say(
    "Best model: ", best$family, " with ", length(best$variables),
    " risk factor(s), BIC = ", round(best$bic, 1)
  )
  if (length(best$variables) == 0) {
    warning("no risk factor met the selection criteria; returning a null model")
  }

  # 5. predictions, relative risk scores, and relative risk values
  grid$outcome_count <- outcome_count
  if (is.null(offset_count)) {
    prediction <- as.numeric(stats::predict(best$fit, type = "response"))
    grid$prediction <- prediction
    grid$relrisk <- prediction / min(prediction)
  } else {
    # rate multiplier exp(b0 + Xb) for every cell (offset set to 0), then
    # expected outcome counts given each cell's denominator
    multiplier_data <- model_data
    multiplier_data$.log_offset <- 0
    multiplier <- as.numeric(
      stats::predict(best$fit, newdata = multiplier_data, type = "response")
    )
    grid$offset_count <- offset_count
    grid$prediction <- offset_count * multiplier
    grid$relrisk <- multiplier / min(multiplier)
  }

  coefs <- stats::coef(summary(best$fit))
  selected <- meta[match(best$variables, meta$variable), ]
  selected$coefficient <- coefs[best$variables, "Estimate"]
  selected$p_value <- coefs[best$variables, ncol(coefs)]
  selected$rrv <- exp(selected$coefficient)
  selected <- selected[order(-selected$rrv), ]
  rownames(selected) <- NULL

  adjustments <- NULL
  if (length(always) > 0) {
    idx <- match(always, rownames(coefs)) # NA if aliased/dropped by glm
    adjustments <- data.frame(
      covariate = always,
      coefficient = coefs[idx, "Estimate"],
      p_value = coefs[idx, ncol(coefs)],
      row.names = NULL
    )
  }

  structure(
    list(
      best_model = best$fit,
      family = best$family,
      bic = best$bic,
      selected = selected,
      adjustments = adjustments,
      grid = grid,
      intercept = unname(stats::coef(best$fit)["(Intercept)"]),
      data = model_data,
      meta = meta,
      has_offset = !is.null(offset_count),
      call = match.call()
    ),
    class = "rtm"
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.rtm <- function(x, ...) {
  cat("Risk Terrain Model\n")
  cat(
    "Distribution:",
    if (x$family == "nb") "negative binomial" else "Poisson",
    sprintf("(BIC = %.1f)\n", x$bic)
  )
  cat("Cells:", nrow(x$grid), " Events:", sum(x$grid$outcome_count))
  if (isTRUE(x$has_offset)) {
    cat("  Denominator:", sum(x$grid$offset_count))
  }
  cat("\n")
  if (isTRUE(x$has_offset)) {
    cat("Rate model: RRVs are adjusted rate ratios (outcome per denominator)\n")
  }
  cat("\n")

  if (nrow(x$selected) == 0) {
    cat("No significant risk factors were found.\n")
  } else {
    cat("Selected risk factors:\n")
    tab <- x$selected[, c(
      "factor", "type", "operation", "spatial_influence", "rrv", "coefficient", "p_value"
    )]
    tab$rrv <- round(tab$rrv, 3)
    tab$coefficient <- round(tab$coefficient, 4)
    tab$p_value <- signif(tab$p_value, 3)
    print(tab, row.names = FALSE)
    cat("\nIntercept:", round(x$intercept, 4), "\n")
    cat(
      if (isTRUE(x$has_offset)) "Relative rates range from" else "Relative risk scores range from",
      round(min(x$grid$relrisk), 2), "to", round(max(x$grid$relrisk), 2), "\n"
    )
  }

  if (!is.null(x$adjustments)) {
    cat("\nAdjustment covariates (included in every model, no RRV):\n")
    adj <- x$adjustments
    adj$coefficient <- round(adj$coefficient, 4)
    adj$p_value <- signif(adj$p_value, 3)
    print(adj, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.rtm <- function(object, ...) {
  print(object)
  cat("\nModel summary:\n")
  print(summary(object$best_model))
  invisible(object)
}

#' Plot the risk terrain map
#'
#' @param x An `rtm` object.
#' @param what Cell value to map: `"relrisk"` (relative risk score, default),
#'   `"prediction"` (expected event count), or `"outcome_count"`.
#' @param ... Passed on to [plot.sf()][sf::plot_sf].
#' @export
plot.rtm <- function(x, what = c("relrisk", "prediction", "outcome_count"), ...) {
  what <- match.arg(what)
  plot(x$grid[what], border = NA, ...)
}

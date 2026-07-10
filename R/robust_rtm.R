#' Cluster-robust inference for a risk terrain model
#'
#' Recomputes the standard errors, p-values, and confidence intervals of a
#' fitted [rtm()] model with a cluster-robust (sandwich) variance estimator,
#' grouping cells into square spatial clusters. Point estimates (and hence
#' RRVs and the risk map) are unchanged; only the uncertainty statements are
#' re-estimated to tolerate correlation between nearby cells.
#'
#' The model-based standard errors assume independent cells. Under residual
#' spatial autocorrelation (see [moran_rtm()]) they are too small and the
#' reported p-values anti-conservative. The cluster-robust estimator allows
#' arbitrary correlation *within* clusters while assuming independence
#' *between* them, so clusters should be large relative to the residual
#' correlation range; by default the cluster side is twice the largest
#' candidate spatial influence distance, so that cells exposed to the same
#' feature usually share a cluster. Because results can be sensitive to
#' this choice, report the robust intervals for two or three cluster sizes
#' when autocorrelation is a concern.
#'
#' Note that `cluster_size` is unrelated to `block_length` in [rtm()]:
#' `block_length` is the average *city block* length used to build the
#' ladder of candidate spatial influences, while `cluster_size` (like
#' `cell_size`, a side length in CRS units) defines the spatial tiles used
#' only for variance estimation.
#'
#' Two limitations to keep in mind. First, cluster asymptotics rely on the
#' *number of clusters*; a warning is issued below 30, where robust
#' p-values are themselves unreliable. Second, robust variances address
#' dependence, not selection: variables were still chosen by the data, so
#' even robust p-values inherit the post-selection caveat, and fully honest
#' inference additionally requires re-estimating the selected specification
#' on holdout data. For negative binomial models the dispersion parameter
#' is treated as fixed at its estimate.
#'
#' @param x An `rtm` object.
#' @param cluster_size Side length of the square spatial clusters, in CRS
#'   units. Defaults to twice the largest candidate spatial influence
#'   distance in the model.
#' @param cluster Optional custom cluster assignment overriding the square
#'   tiling: a vector with one value per grid cell (e.g. neighbourhood or
#'   census tract identifiers).
#' @param level Confidence level (default 0.95).
#'
#' @return An object of class `rtm_robust`: a list with `selected` (the
#'   model's selected-variables table with columns `rrv`, `se_robust` and
#'   `se_model` (log scale), `p_robust`, `p_model`, and RRV-scale confidence
#'   limits `ci_lower`/`ci_upper`), `adjustments` (same columns on the
#'   coefficient scale, `NULL` if none), `n_clusters`, `cluster_size` (`NA`
#'   if a custom `cluster` was supplied), and `level`.
#' @references Cameron, A. C., & Miller, D. L. (2015). A practitioner's
#'   guide to cluster-robust inference. *Journal of Human Resources*, 50,
#'   317-372.
#' @seealso [moran_rtm()] to diagnose whether robust inference is needed,
#'   and [correlogram_rtm()] to estimate `cluster_size` from the residuals'
#'   spatial correlation range.
#' @export
robust_rtm <- function(x, cluster_size = NULL, cluster = NULL, level = 0.95) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("robust_rtm() requires the sandwich package; install it with install.packages(\"sandwich\")")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1")
  }

  fit_cells <- if (isTRUE(x$has_offset)) x$grid$offset_count > 0 else rep(TRUE, nrow(x$grid))

  if (is.null(cluster)) {
    if (is.null(cluster_size)) {
      influences <- x$meta$spatial_influence
      if (all(is.na(influences))) {
        stop(
          "no spatial influence distances found to derive a default cluster ",
          "size; supply `cluster_size` or `cluster`"
        )
      }
      cluster_size <- 2 * max(influences, na.rm = TRUE)
    }
    if (!is.numeric(cluster_size) || length(cluster_size) != 1 || cluster_size <= 0) {
      stop("`cluster_size` must be a single positive number")
    }
    ctr <- suppressWarnings(
      sf::st_coordinates(sf::st_centroid(sf::st_geometry(x$grid)))
    )
    cluster <- paste(
      floor((ctr[, 1] - min(ctr[, 1])) / cluster_size),
      floor((ctr[, 2] - min(ctr[, 2])) / cluster_size)
    )
  } else {
    if (length(cluster) != nrow(x$grid)) {
      stop("`cluster` must have one value per grid cell (", nrow(x$grid), " cells)")
    }
    if (anyNA(cluster)) {
      stop("`cluster` contains missing values")
    }
    cluster_size <- NA_real_
  }

  cl <- factor(cluster[fit_cells])
  n_clusters <- nlevels(cl)
  if (n_clusters < 2) {
    stop("clustering produced fewer than 2 clusters; decrease `cluster_size`")
  }
  if (n_clusters < 30) {
    warning(
      "only ", n_clusters, " spatial clusters; cluster-robust inference ",
      "is unreliable with few clusters - consider a smaller `cluster_size`"
    )
  }

  vc <- sandwich::vcovCL(x$best_model, cluster = cl)
  est <- stats::coef(x$best_model)
  se_model <- sqrt(diag(stats::vcov(x$best_model)))
  se_robust <- sqrt(diag(vc))
  zq <- stats::qnorm(1 - (1 - level) / 2)

  term_table <- function(terms) {
    idx <- match(terms, names(est))
    data.frame(
      coefficient = est[idx],
      se_model = se_model[idx],
      se_robust = se_robust[idx],
      p_model = 2 * stats::pnorm(-abs(est[idx] / se_model[idx])),
      p_robust = 2 * stats::pnorm(-abs(est[idx] / se_robust[idx])),
      ci_lower = est[idx] - zq * se_robust[idx],
      ci_upper = est[idx] + zq * se_robust[idx],
      row.names = NULL
    )
  }

  selected <- x$selected
  if (nrow(selected) > 0) {
    tab <- term_table(selected$variable)
    selected$se_model <- tab$se_model
    selected$se_robust <- tab$se_robust
    selected$p_model <- tab$p_model
    selected$p_robust <- tab$p_robust
    selected$ci_lower <- exp(tab$ci_lower)
    selected$ci_upper <- exp(tab$ci_upper)
    selected$coefficient <- tab$coefficient
    selected$p_value <- NULL
  }

  adjustments <- NULL
  if (!is.null(x$adjustments)) {
    adjustments <- cbind(
      x$adjustments["covariate"],
      term_table(x$adjustments$covariate)
    )
  }

  structure(
    list(
      selected = selected,
      adjustments = adjustments,
      n_clusters = n_clusters,
      cluster_size = cluster_size,
      level = level
    ),
    class = "rtm_robust"
  )
}

#' @export
print.rtm_robust <- function(x, ...) {
  cat("Cluster-robust inference for risk terrain model\n")
  cat(
    "Clusters:", x$n_clusters,
    if (!is.na(x$cluster_size)) sprintf(" (cluster size %g)", x$cluster_size),
    sprintf(" Confidence level: %g%%\n\n", 100 * x$level)
  )

  if (nrow(x$selected) == 0) {
    cat("No selected risk factors.\n")
  } else {
    cat("Selected risk factors (RRV with robust confidence limits):\n")
    tab <- x$selected[, c(
      "factor", "operation", "spatial_influence",
      "rrv", "ci_lower", "ci_upper", "p_robust", "p_model"
    )]
    tab$rrv <- round(tab$rrv, 3)
    tab$ci_lower <- round(tab$ci_lower, 3)
    tab$ci_upper <- round(tab$ci_upper, 3)
    tab$p_robust <- signif(tab$p_robust, 3)
    tab$p_model <- signif(tab$p_model, 3)
    print(tab, row.names = FALSE)
  }

  if (!is.null(x$adjustments)) {
    cat("\nAdjustment covariates (coefficient scale):\n")
    adj <- x$adjustments[, c("covariate", "coefficient", "se_robust", "p_robust")]
    adj$coefficient <- round(adj$coefficient, 4)
    adj$se_robust <- round(adj$se_robust, 4)
    adj$p_robust <- signif(adj$p_robust, 3)
    print(adj, row.names = FALSE)
  }

  cat(
    "\nNote: robust variances address spatial dependence, not model",
    "selection;\np-values remain conditional on the selected specification.\n"
  )
  invisible(x)
}

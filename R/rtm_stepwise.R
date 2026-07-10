fit_count_model <- function(variables, data, outcome, family, offset_col = NULL) {
  rhs <- if (length(variables) == 0) "1" else paste(variables, collapse = " + ")
  if (!is.null(offset_col)) {
    rhs <- paste0(rhs, " + offset(", offset_col, ")")
  }
  f <- stats::as.formula(paste(outcome, "~", rhs))

  tryCatch(
    suppressWarnings(
      if (family == "poisson") {
        stats::glm(f, data = data, family = stats::poisson())
      } else {
        MASS::glm.nb(f, data = data)
      }
    ),
    error = function(e) NULL
  )
}

# RTMDx validity checks: every risk-factor coefficient must adjust risk in the
# theoretically congruent direction (positive, given the 0/1 and 0/-1 coding)
# and be statistically significant; adjustment covariates in `exempt` may
# take either sign and are not required to be significant
model_is_valid <- function(fit, alpha_level = 0.05, exempt = character(0)) {
  if (is.null(fit) || !fit$converged) {
    return(FALSE)
  }
  coefs <- stats::coef(summary(fit))
  coefs <- coefs[!(rownames(coefs) %in% c("(Intercept)", exempt)), , drop = FALSE]
  if (nrow(coefs) == 0) {
    return(TRUE)
  }
  all(coefs[, "Estimate"] > 0) && all(coefs[, ncol(coefs)] <= alpha_level)
}

#' Bidirectional stepwise selection of risk factors by BIC
#'
#' Builds a parsimonious risk terrain model as in the RTMDx Utility: starting
#' from a null model, variables are added and removed one step at a time,
#' keeping the step with the best (lowest) Bayesian information criterion,
#' until no valid step improves the score. RTMDx's entry rules are enforced:
#' at most one variable (i.e. one operationalization and spatial influence)
#' per risk factor may be in the model, all coefficients must be positive
#' (congruent with the coded aggravating/protective direction), and all
#' variables in a candidate model must be significant at `alpha_level`.
#' Note that, as in RTMDx, significance is a sanity check only; steps are
#' chosen by BIC, so reported p-values should not be read as regular p-values.
#'
#' @param data A data.frame with the outcome counts and candidate binary
#'   variables per grid cell.
#' @param outcome Name of the outcome count column in `data`.
#' @param meta Variable metadata data.frame (columns `variable` and `factor`
#'   at minimum), as attached by [operationalize()]; only variables listed
#'   here are candidates.
#' @param family `"poisson"` or `"nb"` (negative binomial, via
#'   [MASS::glm.nb()]).
#' @param alpha_level Significance level for the validity check (default 0.05).
#' @param offset_col Optional name of a column in `data` holding the model
#'   offset on the log scale (e.g. log of denominator event counts).
#' @param always Character vector of adjustment covariate columns included
#'   in every candidate model (starting with the null model). They do not
#'   compete for selection, may take either coefficient sign, and are not
#'   subject to the significance check.
#' @param verbose Print progress messages (default `FALSE`).
#'
#' @return A list with elements `fit` (the final model object), `variables`
#'   (selected variable names, excluding `always` terms), `bic`, and
#'   `family`, or `NULL` if even the null model could not be fitted.
#' @export
rtm_stepwise <- function(data, outcome, meta,
                         family = c("poisson", "nb"),
                         alpha_level = 0.05,
                         offset_col = NULL,
                         always = character(0),
                         verbose = FALSE) {
  family <- match.arg(family)
  candidates <- meta$variable

  current <- character(0)
  current_fit <- fit_count_model(c(always, current), data, outcome, family, offset_col)
  if (is.null(current_fit)) {
    return(NULL)
  }
  current_bic <- stats::BIC(current_fit)

  repeat {
    # additions: candidate variables whose risk factor is not yet represented
    factors_in <- meta$factor[meta$variable %in% current]
    additions <- candidates[
      !(candidates %in% current) & !(meta$factor %in% factors_in)
    ]
    steps <- c(
      lapply(additions, function(v) c(current, v)),
      lapply(current, function(v) setdiff(current, v))
    )

    best_step <- NULL
    best_fit <- NULL
    best_bic <- current_bic

    for (step in steps) {
      fit <- fit_count_model(c(always, step), data, outcome, family, offset_col)
      if (is.null(fit) || !model_is_valid(fit, alpha_level, exempt = always)) {
        next
      }
      bic <- stats::BIC(fit)
      if (bic < best_bic) {
        best_step <- step
        best_fit <- fit
        best_bic <- bic
      }
    }

    if (is.null(best_step)) {
      break
    }

    if (verbose) {
      changed <- c(setdiff(best_step, current), setdiff(current, best_step))
      action <- if (length(best_step) > length(current)) "+" else "-"
      message(sprintf(
        "  step: %s %s (BIC %.1f -> %.1f)", action, changed, current_bic, best_bic
      ))
    }

    current <- best_step
    current_fit <- best_fit
    current_bic <- best_bic
  }

  list(fit = current_fit, variables = current, bic = current_bic, family = family)
}

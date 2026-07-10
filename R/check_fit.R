#' Goodness-of-fit checks for a risk terrain model
#'
#' Assesses how well the final model describes the observed cell counts,
#' from two complementary angles:
#'
#' * **Pseudo-R-squared.** There is no unique R-squared for count models;
#'   three established analogues are reported, all comparing the fitted
#'   model against an intercept-only null (with the same offset).
#'   `deviance` is the proportional reduction in deviance (Cameron &
#'   Windmeijer, 1996), the closest analogue of the OLS R-squared and the
#'   recommended headline number; `mcfadden` is the log-likelihood
#'   version; `nagelkerke` is included for comparability with published
#'   RTM studies, which report it. Expect *small* values on cell-level
#'   count data even for genuinely good models: with sparse counts, most
#'   cell-to-cell variation is irreducible Poisson noise that no
#'   covariate can explain. The measures are meaningful for comparing
#'   models on the *same* data, not against OLS conventions; for
#'   predictive performance use [validate_rtm()].
#' * **Posterior-predictive checks.** `nsim` replicate datasets are
#'   simulated from the fitted model (Poisson, or negative binomial with
#'   the estimated dispersion) and summary statistics of the observed
#'   counts are compared against their simulated distributions: the share
#'   of zero cells (are the zeros ordinary sampling zeros? see the
#'   `exclude` argument of [rtm()] for structural zeros), the dispersion
#'   index (variance/mean; residual overdispersion), and the maximum cell
#'   count (extreme hot cells). A small two-sided p-value means the model
#'   cannot reproduce that feature of the data.
#'
#' For offset models all quantities use the fitted (non-zero denominator)
#' cells. If the model includes adjustment covariates, they count towards
#' the "explained" side of the pseudo-R-squared, which then reflects the
#' full specification rather than the selected risk factors alone.
#'
#' @param x An `rtm` object.
#' @param nsim Number of simulated replicate datasets (default 500).
#' @param level Level for the reported posterior-predictive intervals
#'   (default 0.95).
#'
#' @return An object of class `rtm_fit_check`: a list with `pseudo_r2`
#'   (named numeric vector: `deviance`, `mcfadden`, `nagelkerke`),
#'   `checks` (data.frame with the observed statistic, the simulated
#'   interval bounds `sim_lo`/`sim_hi`, the simulated median, and the
#'   two-sided permutation p-value for `zero_share`, `dispersion`, and
#'   `max_count`), `counts`/`mu`/`theta` (observed counts, fitted means,
#'   and NB dispersion, used by [plot.rtm_fit_check()]), `family`,
#'   `n_cells`, `nsim`, and `level`.
#' @references Cameron, A. C., & Windmeijer, F. A. G. (1996). R-squared
#'   measures for count data regression models with applications to
#'   health-care utilization. *Journal of Business & Economic
#'   Statistics*, 14, 209-220.
#'
#'   Nagelkerke, N. J. D. (1991). A note on a general definition of the
#'   coefficient of determination. *Biometrika*, 78, 691-692.
#' @seealso [validate_rtm()] for out-of-sample predictive accuracy,
#'   [moran_rtm()] for residual spatial structure.
#' @export
check_fit <- function(x, nsim = 500, level = 0.95) {
  if (!inherits(x, "rtm")) {
    stop("`x` must be an rtm object")
  }
  if (nsim < 2) {
    stop("`nsim` must be at least 2")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1")
  }

  fit <- x$best_model
  fit_cells <- if (isTRUE(x$has_offset)) x$grid$offset_count > 0 else rep(TRUE, nrow(x$grid))
  y <- x$data$outcome_count[fit_cells]
  mu <- as.numeric(stats::fitted(fit))
  if (length(mu) != length(y)) {
    stop(
      "fitted values (", length(mu), ") do not match the fitted cells (",
      length(y), "); was the object modified?"
    )
  }
  n <- length(y)

  # pseudo-R-squared against the intercept-only null with the same offset;
  # the object's null deviance uses the model's own dispersion, the null
  # log-likelihood comes from refitting the null model
  offset_col <- if (isTRUE(x$has_offset)) ".log_offset"
  null_fit <- fit_count_model(
    character(0), x$data[fit_cells, , drop = FALSE],
    outcome = "outcome_count", family = x$family, offset_col = offset_col
  )
  r2_deviance <- 1 - stats::deviance(fit) / fit$null.deviance
  r2_mcfadden <- NA_real_
  r2_nagelkerke <- NA_real_
  if (!is.null(null_fit)) {
    ll1 <- as.numeric(stats::logLik(fit))
    ll0 <- as.numeric(stats::logLik(null_fit))
    r2_mcfadden <- 1 - ll1 / ll0
    cox_snell <- 1 - exp(2 * (ll0 - ll1) / n)
    r2_nagelkerke <- cox_snell / (1 - exp(2 * ll0 / n))
  } else {
    warning("the null model could not be fitted; likelihood-based pseudo-R2 values are NA")
  }

  # posterior-predictive checks: simulate replicate counts from the fitted
  # model and compare summary statistics
  theta <- if (x$family == "nb") fit$theta
  simulate_counts <- if (x$family == "nb") {
    function() stats::rnbinom(n, mu = mu, size = theta)
  } else {
    function() stats::rpois(n, mu)
  }
  stat_fun <- function(v) {
    c(
      zero_share = mean(v == 0),
      dispersion = stats::var(v) / mean(v),
      max_count = max(v)
    )
  }

  observed <- stat_fun(y)
  sims <- vapply(seq_len(nsim), function(s) stat_fun(simulate_counts()), observed)

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  checks <- do.call(rbind, lapply(names(observed), function(nm) {
    s <- sims[nm, ]
    q <- stats::quantile(s, probs, names = FALSE)
    p_hi <- (1 + sum(s >= observed[nm])) / (nsim + 1)
    p_lo <- (1 + sum(s <= observed[nm])) / (nsim + 1)
    data.frame(
      statistic = nm,
      observed = unname(observed[nm]),
      sim_lo = q[1],
      sim_median = q[2],
      sim_hi = q[3],
      p_value = min(1, 2 * min(p_hi, p_lo))
    )
  }))

  structure(
    list(
      pseudo_r2 = c(
        deviance = r2_deviance,
        mcfadden = r2_mcfadden,
        nagelkerke = r2_nagelkerke
      ),
      checks = checks,
      counts = y,
      mu = mu,
      theta = theta,
      family = x$family,
      n_cells = n,
      nsim = nsim,
      level = level
    ),
    class = "rtm_fit_check"
  )
}

# expected frequency of each count value under the fitted model, summing
# each cell's Poisson/NB probability at its own fitted mean
expected_frequencies <- function(mu, theta, k) {
  vapply(k, function(kk) {
    if (is.null(theta)) {
      sum(stats::dpois(kk, mu))
    } else {
      sum(stats::dnbinom(kk, mu = mu, size = theta))
    }
  }, numeric(1))
}

# randomized quantile residuals (Dunn & Smyth, 1996): uniform draws between
# the fitted CDF just below and at the observed count, mapped to normal
quantile_residuals <- function(y, mu, theta = NULL) {
  cdf <- function(q) {
    if (is.null(theta)) {
      stats::ppois(q, mu)
    } else {
      stats::pnbinom(q, mu = mu, size = theta)
    }
  }
  u <- stats::runif(length(y), cdf(y - 1), cdf(y))
  stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}

#' Plot goodness-of-fit diagnostics for a risk terrain model
#'
#' Draws the two standard visual fit checks for count models from a
#' [check_fit()] result:
#'
#' * **Hanging rootogram** (Kleiber & Zeileis, 2016): for each count value
#'   0, 1, 2, ... the square root of the observed cell frequency hangs as
#'   a bar from the square root of the model-expected frequency. Bars
#'   dangling below the zero line are counts the model *under*-predicts;
#'   bars stopping short of it are counts it *over*-predicts. Zero excess
#'   shows as the 0-bar punching below the line; overdispersion as a
#'   sagging middle with heavy tails.
#' * **Quantile-residual Q-Q plot** (Dunn & Smyth, 1996): randomized
#'   quantile residuals are standard normal when the model is correct, so
#'   departures from the diagonal read like an ordinary Q-Q plot. The
#'   randomization draws from the RNG; call `set.seed()` first for a
#'   reproducible figure.
#'
#' @param x An `rtm_fit_check` object from [check_fit()].
#' @param which Plots to draw: `"rootogram"`, `"qq"`, or both (default;
#'   drawn side by side).
#' @param max_count Largest count value shown in the rootogram (default:
#'   the maximum observed count).
#' @param ... Passed on to [plot()].
#' @references Kleiber, C., & Zeileis, A. (2016). Visualizing count data
#'   regressions using rootograms. *The American Statistician*, 70,
#'   296-303.
#'
#'   Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
#'   *Journal of Computational and Graphical Statistics*, 5, 236-244.
#' @export
plot.rtm_fit_check <- function(x, which = c("rootogram", "qq"),
                               max_count = NULL, ...) {
  which <- match.arg(which, several.ok = TRUE)
  if (length(which) > 1) {
    op <- graphics::par(mfrow = c(1, length(which)))
    on.exit(graphics::par(op))
  }

  if ("rootogram" %in% which) {
    k <- 0:(max_count %||% max(x$counts))
    obs <- vapply(k, function(kk) sum(x$counts == kk), numeric(1))
    expd <- expected_frequencies(x$mu, x$theta, k)
    top <- sqrt(expd)
    bottom <- top - sqrt(obs)
    plot(
      k, top,
      type = "n", xlab = "Count per cell", ylab = "sqrt(frequency)",
      ylim = range(0, top, bottom), main = "Hanging rootogram", ...
    )
    graphics::rect(k - 0.4, bottom, k + 0.4, top, col = "grey85", border = "grey40")
    graphics::lines(k, top, type = "o", pch = 16, col = "red3")
    graphics::abline(h = 0, lty = 2)
    graphics::mtext("bar below 0 = under-predicted count", side = 3, cex = 0.7)
  }

  if ("qq" %in% which) {
    r <- sort(quantile_residuals(x$counts, x$mu, x$theta))
    q <- stats::qnorm(stats::ppoints(length(r)))
    plot(
      q, r,
      xlab = "Theoretical normal quantiles", ylab = "Quantile residuals",
      main = "Quantile-residual Q-Q", pch = 16,
      col = grDevices::adjustcolor("black", alpha.f = 0.4), ...
    )
    graphics::abline(0, 1, col = "red3")
  }

  invisible(x)
}

#' @export
print.rtm_fit_check <- function(x, ...) {
  cat("Goodness of fit of the risk terrain model\n")
  cat(
    "Family:", if (x$family == "nb") "negative binomial" else "Poisson",
    " Cells:", x$n_cells, " Simulations:", x$nsim, "\n\n"
  )

  cat("Pseudo-R-squared (vs intercept-only null):\n")
  r2 <- round(x$pseudo_r2, 3)
  cat(sprintf(
    "  deviance (Cameron-Windmeijer): %s\n  McFadden: %s\n  Nagelkerke: %s\n",
    r2["deviance"], r2["mcfadden"], r2["nagelkerke"]
  ))
  cat(
    "  Note: small values are expected for sparse cell counts (most\n",
    "  cell-level variation is irreducible count noise); compare models\n",
    "  on the same data rather than against OLS conventions.\n",
    sep = ""
  )

  cat("\nPosterior-predictive checks (observed vs simulated from the model):\n")
  tab <- x$checks
  tab$observed <- signif(tab$observed, 4)
  tab$sim_lo <- signif(tab$sim_lo, 4)
  tab$sim_median <- signif(tab$sim_median, 4)
  tab$sim_hi <- signif(tab$sim_hi, 4)
  tab$p_value <- signif(tab$p_value, 3)
  print(tab, row.names = FALSE)

  flagged <- x$checks$statistic[x$checks$p_value <= 0.05]
  if (length(flagged) > 0) {
    cat("\nThe model cannot reproduce: ", paste(flagged, collapse = ", "), ".\n", sep = "")
    if ("zero_share" %in% flagged) {
      cat(
        "Excess zeros suggest unmasked structural zeros (see `exclude` in\n",
        "?rtm) or a zero-inflated process.\n",
        sep = ""
      )
    }
    if ("dispersion" %in% flagged && x$family == "poisson") {
      cat("Residual overdispersion: the negative binomial model may fit better.\n")
    }
  }
  invisible(x)
}

#' Assign grid cells to cross-validation folds stratified on the outcome
#'
#' Replicates the RTMDx fold-building procedure: cells are sorted by outcome
#' count (high to low) with random tie-breaking, then allocated to folds in
#' consecutive blocks of `nfolds` cells with random order within each block,
#' so that outcome counts are balanced across folds.
#'
#' @param y Integer vector of outcome counts per cell.
#' @param nfolds Number of folds (default 5, as in RTMDx).
#'
#' @return An integer vector of fold ids (1 to `nfolds`), one per cell.
#' @export
make_folds <- function(y, nfolds = 5) {
  n <- length(y)
  ord <- order(-y, stats::runif(n))

  allocation <- as.vector(
    replicate(ceiling(n / nfolds), sample.int(nfolds))
  )[seq_len(n)]

  folds <- integer(n)
  folds[ord] <- allocation
  folds
}

#' Cull candidate variables with penalized Poisson regression
#'
#' Reduces the large set of candidate spatial-influence variables to those
#' worth testing further, guarding against spurious correlations from
#' multiple comparisons. Following RTMDx, an elastic-net penalized Poisson
#' regression is fitted with the L1 (lasso) penalty optimized by
#' cross-validation over folds stratified on the outcome, and coefficients
#' constrained to be non-negative (variables are coded so that their expected
#' effect is positive). Variables with non-zero coefficients survive.
#'
#' RTMDx fits two small fixed L2 (ridge) penalties and picks the one with the
#' best cross-validated likelihood; here the analogous mixing is controlled
#' by trying several `alpha` values in [glmnet::cv.glmnet()] (alpha < 1 adds
#' an L2 component, which helps collinear variables survive together).
#'
#' Fold assignment is randomized (as in RTMDx, where it changes on every
#' invocation), so which variables survive can vary from run to run; call
#' `set.seed()` beforehand for exact reproducibility. Setting `repeats`
#' above 1 additionally stabilizes the result by *repeated*
#' cross-validation: the cross-validated deviance curve is averaged over
#' `repeats` independent fold assignments before the penalty is chosen
#' (the lambda path is computed from the full data, so it is common across
#' repeats), making the culling far less sensitive to any single fold
#' draw. `repeats = 1` reproduces the RTMDx single-assignment behaviour.
#'
#' @param x A data.frame or matrix of candidate binary variables (columns)
#'   per grid cell (rows), e.g. built with [operationalize()].
#' @param y Integer vector of outcome counts per cell.
#' @param alpha Elastic-net mixing values to try; the fit with the best
#'   cross-validated deviance is used.
#' @param nfolds Number of cross-validation folds (default 5).
#' @param offset Optional numeric vector with the model offset on the log
#'   scale (e.g. `log()` of denominator event counts), one value per cell.
#' @param adjust Optional data.frame or matrix of adjustment covariates that
#'   are included unpenalized and sign-unconstrained; they do not compete
#'   for selection and are not returned.
#' @param repeats Number of independent fold assignments to average the
#'   cross-validation over (default 1, as in RTMDx; 5-10 recommended for
#'   stable results).
#'
#' @return A character vector with the names of the surviving candidate
#'   variables.
#' @export
cull_variables <- function(x, y, alpha = c(0.8, 0.95), nfolds = 5, offset = NULL,
                           adjust = NULL, repeats = 1) {
  if (!is.numeric(repeats) || length(repeats) != 1 || repeats < 1) {
    stop("`repeats` must be a single positive integer")
  }
  repeats <- as.integer(repeats)
  x <- as.matrix(x)
  candidates <- colnames(x)

  if (!is.null(adjust)) {
    adjust <- as.matrix(adjust)
    penalty <- c(rep(0, ncol(adjust)), rep(1, ncol(x)))
    lower <- c(rep(-Inf, ncol(adjust)), rep(0, ncol(x)))
    x <- cbind(adjust, x)
  } else {
    penalty <- rep(1, ncol(x))
    lower <- rep(0, ncol(x))
  }

  # one fold assignment per repeat, shared across the alpha values so that
  # alphas are compared on equal folds within each repeat
  fits <- lapply(seq_len(repeats), function(r) {
    foldid <- make_folds(y, nfolds)
    lapply(alpha, function(a) {
      tryCatch(
        glmnet::cv.glmnet(x, y,
          family = "poisson", alpha = a,
          lower.limits = lower, penalty.factor = penalty,
          foldid = foldid, offset = offset
        ),
        error = function(e) NULL
      )
    })
  })

  # per alpha, average the CV deviance curves over the repeats; the lambda
  # path is data-determined and hence shared, up to early stopping, so the
  # curves are aligned on their common leading lambdas
  best_cvm <- Inf
  best_fit <- NULL
  best_lambda <- NULL
  for (j in seq_along(alpha)) {
    afits <- Filter(Negate(is.null), lapply(fits, `[[`, j))
    if (length(afits) == 0) next
    k <- min(vapply(afits, function(f) length(f$lambda), integer(1)))
    avg <- rowMeans(
      matrix(vapply(afits, function(f) f$cvm[seq_len(k)], numeric(k)), nrow = k)
    )
    i <- which.min(avg)
    if (avg[i] < best_cvm) {
      best_cvm <- avg[i]
      best_fit <- afits[[1]]
      best_lambda <- afits[[1]]$lambda[i]
    }
  }

  if (is.null(best_fit)) {
    warning("penalized regression failed to converge; keeping all variables")
    return(candidates)
  }

  coefs <- stats::coef(best_fit, s = best_lambda)
  coefs <- coefs[rownames(coefs) %in% candidates, 1]
  names(coefs)[coefs != 0]
}

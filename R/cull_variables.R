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
#'
#' @return A character vector with the names of the surviving candidate
#'   variables.
#' @export
cull_variables <- function(x, y, alpha = c(0.8, 0.95), nfolds = 5, offset = NULL,
                           adjust = NULL) {
  x <- as.matrix(x)
  candidates <- colnames(x)
  foldid <- make_folds(y, nfolds)

  if (!is.null(adjust)) {
    adjust <- as.matrix(adjust)
    penalty <- c(rep(0, ncol(adjust)), rep(1, ncol(x)))
    lower <- c(rep(-Inf, ncol(adjust)), rep(0, ncol(x)))
    x <- cbind(adjust, x)
  } else {
    penalty <- rep(1, ncol(x))
    lower <- rep(0, ncol(x))
  }

  fits <- lapply(alpha, function(a) {
    tryCatch(
      glmnet::cv.glmnet(x, y,
        family = "poisson", alpha = a,
        lower.limits = lower, penalty.factor = penalty,
        foldid = foldid, offset = offset
      ),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)

  if (length(fits) == 0) {
    warning("penalized regression failed to converge; keeping all variables")
    return(candidates)
  }

  cv_dev <- vapply(fits, function(f) min(f$cvm), numeric(1))
  best <- fits[[which.min(cv_dev)]]

  coefs <- stats::coef(best, s = "lambda.min")
  coefs <- coefs[rownames(coefs) %in% candidates, 1]
  names(coefs)[coefs != 0]
}

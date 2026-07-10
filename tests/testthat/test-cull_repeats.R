# candidate matrix with one true predictor and noise columns
cull_test_data <- function(n = 400, seed = 5) {
  set.seed(seed)
  x <- data.frame(
    real = rbinom(n, 1, 0.3),
    noise1 = rbinom(n, 1, 0.3),
    noise2 = rbinom(n, 1, 0.3),
    noise3 = rbinom(n, 1, 0.3)
  )
  y <- rpois(n, exp(-1 + 1.5 * x$real))
  list(x = x, y = y)
}

test_that("repeats = 1 reproduces the single-assignment culling exactly", {
  d <- cull_test_data()

  # old behaviour: one fold assignment shared across alphas, lambda.min of
  # the best alpha by min cvm
  set.seed(1)
  foldid <- make_folds(d$y, 5)
  fits <- lapply(c(0.8, 0.95), function(a) {
    glmnet::cv.glmnet(
      as.matrix(d$x), d$y,
      family = "poisson", alpha = a,
      lower.limits = rep(0, ncol(d$x)),
      penalty.factor = rep(1, ncol(d$x)), foldid = foldid
    )
  })
  best <- fits[[which.min(vapply(fits, function(f) min(f$cvm), numeric(1)))]]
  coefs <- stats::coef(best, s = "lambda.min")
  old_kept <- rownames(coefs)[-1][coefs[-1, 1] != 0]

  set.seed(1)
  new_kept <- cull_variables(d$x, d$y, repeats = 1)
  expect_equal(sort(new_kept), sort(old_kept))
})

test_that("repeated culling keeps the true predictor and is deterministic", {
  d <- cull_test_data()

  set.seed(2)
  kept_a <- cull_variables(d$x, d$y, repeats = 5)
  set.seed(2)
  kept_b <- cull_variables(d$x, d$y, repeats = 5)

  expect_identical(kept_a, kept_b)
  expect_true("real" %in% kept_a)
  expect_true(all(kept_a %in% names(d$x)))
})

test_that("repeated culling is less fold-sensitive than a single assignment", {
  d <- cull_test_data()

  runs <- function(repeats) {
    vapply(1:12, function(s) {
      set.seed(100 + s)
      paste(sort(cull_variables(d$x, d$y, repeats = repeats)), collapse = "+")
    }, character(1))
  }

  # averaging over fold draws cannot increase the number of distinct
  # culling outcomes across seeds
  expect_lte(length(unique(runs(8))), length(unique(runs(1))))
})

test_that("cull_variables validates repeats and rtm passes cull_repeats through", {
  d <- cull_test_data()
  expect_error(cull_variables(d$x, d$y, repeats = 0), "repeats")

  s <- synthetic_rtm_data(seed = 42)
  set.seed(1)
  fit <- rtm(
    outcome = s$crimes,
    factors = list(bars = s$bars, parks = s$parks),
    boundary = s$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    cull_repeats = 3,
    verbose = FALSE
  )
  expect_s3_class(fit, "rtm")
  expect_true("bars" %in% fit$selected$factor)
})

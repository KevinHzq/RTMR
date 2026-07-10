test_that("make_folds balances outcome counts across folds", {
  y <- c(rep(0, 95), rep(5, 5))
  folds <- make_folds(y, nfolds = 5)
  expect_equal(sort(unique(folds)), 1:5)
  # the five high-count cells land in five different folds
  expect_equal(sort(folds[y == 5]), 1:5)
})

test_that("rtm recovers a planted risk factor end to end", {
  d <- synthetic_rtm_data()

  set.seed(1)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars, parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    max_blocks = 3,
    verbose = FALSE
  )

  expect_s3_class(fit, "rtm")
  expect_true(fit$family %in% c("poisson", "nb"))
  expect_true("bars" %in% fit$selected$factor)

  # at most one variable per risk factor
  expect_false(any(duplicated(fit$selected$factor)))

  # all coefficients positive and significant (RTMDx validity rules)
  expect_true(all(fit$selected$coefficient > 0))
  expect_true(all(fit$selected$p_value <= 0.05))
  expect_true(all(fit$selected$rrv > 1))

  # relative risk scores start at 1
  expect_equal(min(fit$grid$relrisk), 1)
  expect_true(all(c("outcome_count", "prediction", "relrisk") %in% names(fit$grid)))
  expect_equal(sum(fit$grid$outcome_count), length(d$crimes))

  expect_output(print(fit), "Risk Terrain Model")
  expect_output(print(fit), "bars")
})

test_that("rtm accepts per-factor specifications", {
  d <- synthetic_rtm_data()

  set.seed(2)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(
      bars = list(data = d$bars, operation = "both", increment = "half"),
      parks = list(data = d$parks, max_blocks = 2, type = "protective")
    ),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    cull = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "rtm")
  bars_meta <- fit$meta[fit$meta$factor == "bars", ]
  parks_meta <- fit$meta[fit$meta$factor == "parks", ]
  expect_true(all(c("proximity", "density") %in% bars_meta$operation))
  expect_true(all(parks_meta$spatial_influence <= 600))
  expect_equal(unique(parks_meta$type), "protective")
})

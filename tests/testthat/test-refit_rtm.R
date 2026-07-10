fit_synthetic <- function(events_per_bar = 8) {
  d <- synthetic_rtm_data(seed = 42, events_per_bar = events_per_bar)
  set.seed(1)
  list(
    fit = rtm(
      outcome = d$crimes,
      factors = list(bars = d$bars, parks = d$parks),
      boundary = d$boundary,
      cell_size = 150,
      block_length = 300,
      operation = "proximity",
      verbose = FALSE
    ),
    d = d
  )
}

test_that("split_events partitions the layer reproducibly", {
  d <- synthetic_rtm_data(seed = 42)
  set.seed(7)
  sp <- split_events(d$crimes)
  n <- length(d$crimes)

  expect_equal(length(sp$train) + length(sp$test), n)
  # the split is a partition: recombining recovers every original point
  both <- c(sf::st_geometry(sp$train), sf::st_geometry(sp$test))
  expect_setequal(
    apply(sf::st_coordinates(both), 1, paste, collapse = "_"),
    apply(sf::st_coordinates(d$crimes), 1, paste, collapse = "_")
  )
  # roughly balanced at prop = 0.5
  expect_gt(length(sp$train), n * 0.25)
  expect_lt(length(sp$train), n * 0.75)

  # sf input returns sf halves with rows preserved
  crimes_sf <- sf::st_sf(id = seq_len(n), geometry = d$crimes)
  set.seed(7)
  sp2 <- split_events(crimes_sf)
  expect_s3_class(sp2$train, "sf")
  expect_setequal(c(sp2$train$id, sp2$test$id), seq_len(n))

  expect_error(split_events(d$crimes, prop = 1), "strictly between")
  expect_error(split_events(1:10), "sf or sfc")
})

test_that("refitting on the same events reproduces the original model", {
  s <- fit_synthetic()
  refit <- refit_rtm(s$fit, s$d$crimes)

  expect_s3_class(refit, "rtm_refit")
  expect_equal(refit$selected$rrv_refit, refit$selected$rrv_original, tolerance = 1e-8)
  expect_equal(refit$intercept, s$fit$intercept, tolerance = 1e-8)
  expect_equal(refit$family, s$fit$family)
  expect_true(refit$converged)
  expect_output(print(refit), "original vs refitted")
})

test_that("thinned-split refit recovers the effect without selection", {
  # more events so both halves support estimation
  s <- fit_synthetic(events_per_bar = 30)
  set.seed(8)
  sp <- split_events(s$d$crimes)

  set.seed(9)
  fit_train <- rtm(
    outcome = sp$train,
    factors = list(bars = s$d$bars, parks = s$d$parks),
    boundary = s$d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
  skip_if(nrow(fit_train$selected) == 0, "training half selected no factors")

  refit <- refit_rtm(fit_train, sp$test)
  # bars is a genuine strong effect: it should replicate on the test half
  bars_row <- refit$selected[refit$selected$factor == "bars", ]
  expect_equal(nrow(bars_row), 1)
  expect_gt(bars_row$rrv_refit, 1)
  expect_lte(bars_row$p_refit, 0.05)
  # CI brackets its own point estimate
  expect_true(all(refit$selected$ci_lower <= refit$selected$rrv_refit))
  expect_true(all(refit$selected$ci_upper >= refit$selected$rrv_refit))
})

test_that("refit_rtm handles offset models, reusing or replacing the denominator", {
  d <- synthetic_offset_data(seed = 7)
  set.seed(5)
  fit <- rtm(
    outcome = d$fatal_events,
    offset = d$all_events,
    factors = list(bars = d$bars, naloxone = list(data = d$naloxone, type = "protective")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )

  # same-period style: reuse the original denominator
  same <- refit_rtm(fit, d$fatal_events)
  expect_equal(same$selected$rrv_refit, same$selected$rrv_original, tolerance = 1e-8)

  # temporal holdout style: period-2 numerator and denominator
  hold <- suppressWarnings(refit_rtm(fit, d$fatal_events2, new_offset = d$all_events2))
  expect_true(all(is.finite(hold$selected$rrv_refit)))
  expect_true(all(hold$selected$p_refit > 0 & hold$selected$p_refit <= 1))
  # the generating process is stationary, so the naloxone effect replicates
  nal <- hold$selected[hold$selected$factor == "naloxone", ]
  if (nrow(nal) == 1) {
    expect_lte(nal$p_refit, 0.05)
  }
})

test_that("refit_rtm validates its inputs", {
  s <- fit_synthetic()
  far_away <- sf::st_sfc(sf::st_point(c(1e6, 1e6)), crs = 32610)
  expect_error(refit_rtm(s$fit, far_away), "no new events")
  expect_error(refit_rtm(s$fit, s$d$crimes, new_offset = rep(1, 400)), "only applies")
  expect_error(refit_rtm(s$fit, s$d$crimes, level = 0), "level")
  expect_error(refit_rtm(list(), s$d$crimes), "rtm object")
})

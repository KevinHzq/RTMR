# fit on period 1, validate on period 2 events from the same process
fit_and_holdout <- function() {
  d <- synthetic_rtm_data(seed = 42)
  set.seed(1)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars, parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )

  # second period: new events around the same bars (same generating process)
  set.seed(99)
  bars <- sf::st_coordinates(d$bars)
  new_events <- do.call(rbind, lapply(seq_len(nrow(bars)), function(i) {
    n <- stats::rpois(1, 8)
    cbind(stats::rnorm(n, bars[i, 1], 100), stats::rnorm(n, bars[i, 2], 100))
  }))
  new_events <- new_events[
    new_events[, 1] >= 0 & new_events[, 1] <= 3000 &
      new_events[, 2] >= 0 & new_events[, 2] <= 3000, ,
    drop = FALSE
  ]
  new_events <- sf::st_sfc(
    lapply(seq_len(nrow(new_events)), function(i) sf::st_point(new_events[i, ])),
    crs = 32610
  )

  list(fit = fit, new_events = new_events)
}

test_that("validate_rtm reports capture, PAI, and AUC on holdout events", {
  v <- fit_and_holdout()
  val <- validate_rtm(v$fit, v$new_events)

  expect_s3_class(val, "rtm_validation")
  expect_equal(val$n_events, length(v$new_events))
  expect_equal(sum(val$counts), val$n_events)

  # the model was fit on the true generating process, so it should
  # discriminate well and capture events far above the area-wide rate
  expect_gt(val$auc, 0.7)
  expect_gt(val$capture$pai[val$capture$area_pct == 10], 2)

  # capture is monotone in flagged area; PAI cannot exceed 100 / area_pct
  expect_true(all(diff(val$capture$pct_events) >= 0))
  expect_true(all(val$capture$pct_events <= 100))
  expect_true(all(val$capture$pai <= 100 / val$capture$area_pct))

  # ROC curve runs from (0, 0) to (1, 1) and is non-decreasing
  expect_equal(val$roc$fpr[1], 0)
  expect_equal(val$roc$tpr[nrow(val$roc)], 1)
  expect_true(all(diff(val$roc$tpr) >= 0))
  expect_true(all(diff(val$roc$fpr) >= 0))

  expect_output(print(val), "AUC")
})

test_that("validate_rtm AUC matches a brute-force pairwise computation", {
  v <- fit_and_holdout()
  val <- validate_rtm(v$fit, v$new_events)

  score <- v$fit$grid$relrisk
  pos <- score[val$counts > 0]
  neg <- score[val$counts == 0]
  pairs <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  expect_equal(val$auc, mean(pairs), tolerance = 1e-10)
})

test_that("validate_rtm errors when no events fall in the grid", {
  v <- fit_and_holdout()
  far_away <- sf::st_sfc(sf::st_point(c(1e6, 1e6)), crs = 32610)
  expect_error(validate_rtm(v$fit, far_away), "no holdout events")
})

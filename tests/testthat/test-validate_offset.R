fit_offset_model <- function(d) {
  set.seed(3)
  rtm(
    outcome = d$fatal_events,
    offset = d$all_events,
    factors = list(
      bars = d$bars,
      naloxone = list(data = d$naloxone, type = "protective")
    ),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
}

test_that("denominator-aware validation reports adjusted PAI and calibration", {
  d <- synthetic_offset_data()
  fit <- fit_offset_model(d)

  val <- validate_rtm(fit, d$fatal_events2, new_offset = d$all_events2)

  expect_s3_class(val, "rtm_validation")
  expect_true(val$denominator_aware)
  expect_true(all(c("n_denom", "pct_denom", "adj_pai") %in% names(val$capture)))

  # the lethality model ranks by relative fatality: holdout fatal events
  # should concentrate in flagged cells beyond the denominator share
  expect_gt(val$capture$adj_pai[val$capture$area_pct == 20], 1)

  # denominator shares are monotone and bounded
  expect_true(all(diff(val$capture$pct_denom) >= 0))
  expect_true(all(val$capture$pct_denom <= 100))

  # calibration against the holdout denominator: same generating process,
  # so observed/expected should be near 1
  expect_gt(val$observed_expected, 0.6)
  expect_lt(val$observed_expected, 1.5)

  expect_output(print(val), "adj_pai")
})

test_that("plain validation still works for offset models (area-based)", {
  d <- synthetic_offset_data()
  fit <- fit_offset_model(d)

  val <- validate_rtm(fit, d$fatal_events2)
  expect_false(val$denominator_aware)
  expect_false("adj_pai" %in% names(val$capture))
})

test_that("new_offset input is checked", {
  d <- synthetic_offset_data()
  fit <- fit_offset_model(d)

  expect_error(
    validate_rtm(fit, d$fatal_events2, new_offset = c(1, 2, 3)),
    "one value per grid cell"
  )

  # new_offset is rejected for models without an offset
  d2 <- synthetic_rtm_data()
  set.seed(1)
  fit2 <- rtm(
    outcome = d2$crimes, factors = list(bars = d2$bars),
    boundary = d2$boundary, cell_size = 150, block_length = 300,
    operation = "proximity", verbose = FALSE
  )
  expect_error(
    validate_rtm(fit2, d2$crimes, new_offset = d2$crimes),
    "only applies to models fitted with an offset"
  )
})

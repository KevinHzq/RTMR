test_that("rtm with an event-count offset recovers a lethality factor", {
  d <- synthetic_offset_data()

  set.seed(3)
  # fatal events are a subset of all events, so no cell can have an outcome
  # without a denominator and no exclusion warning should fire
  expect_no_warning(
    fit <- rtm(
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
  )

  expect_s3_class(fit, "rtm")
  expect_true(fit$has_offset)

  # the lethality-driving factor is found and coded protective
  expect_true("naloxone" %in% fit$selected$factor)
  nalox <- fit$selected[fit$selected$factor == "naloxone", ]
  expect_equal(nalox$type, "protective")
  expect_gt(nalox$rrv, 1)

  # grid outputs: denominator column, predictions bounded sensibly, relrisk
  # defined for every cell (including zero-denominator ones) starting at 1
  expect_true(all(c("offset_count", "prediction", "relrisk") %in% names(fit$grid)))
  expect_equal(sum(fit$grid$offset_count), length(d$all_events))
  expect_equal(min(fit$grid$relrisk), 1)
  expect_true(all(fit$grid$prediction[fit$grid$offset_count == 0] == 0))
  expect_true(all(is.finite(fit$grid$relrisk)))

  # expected counts track the denominator: total prediction is close to
  # total observed fatal events (model is calibrated in aggregate)
  expect_equal(
    sum(fit$grid$prediction),
    sum(fit$grid$outcome_count[fit$grid$offset_count > 0]),
    tolerance = 0.15
  )

  expect_output(print(fit), "Rate model")
  expect_output(print(fit), "Denominator")
})

test_that("rtm accepts a numeric per-cell offset", {
  d <- synthetic_offset_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  pop <- count_points(d$all_events, grid) + 1 # strictly positive denominator

  set.seed(4)
  fit <- suppressWarnings(rtm(
    outcome = d$fatal_events,
    offset = pop,
    factors = list(naloxone = list(data = d$naloxone, type = "protective")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  ))

  expect_s3_class(fit, "rtm")
  expect_true(fit$has_offset)
  expect_equal(fit$grid$offset_count, pop)

  # wrong length errors
  expect_error(
    rtm(
      outcome = d$fatal_events, offset = pop[-1],
      factors = list(naloxone = d$naloxone),
      boundary = d$boundary, cell_size = 150, block_length = 300,
      verbose = FALSE
    ),
    "one value per grid cell"
  )
})

test_that("models without an offset are unchanged", {
  d <- synthetic_rtm_data()
  set.seed(1)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
  expect_false(fit$has_offset)
  expect_false("offset_count" %in% names(fit$grid))
  expect_equal(
    fit$grid$prediction,
    as.numeric(predict(fit$best_model, type = "response"))
  )
})

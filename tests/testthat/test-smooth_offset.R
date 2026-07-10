test_that("smooth_offset preserves the total denominator count", {
  d <- synthetic_offset_data()
  grid <- create_grid(d$boundary, cellsize = 150)

  sm <- smooth_offset(d$all_events, grid, bandwidth = 450)

  expect_length(sm, nrow(grid))
  expect_true(all(sm >= 0))
  expect_equal(sum(sm), sum(count_points(d$all_events, grid)))

  # smoothing spreads support: fewer zero-denominator cells than raw counts
  raw <- count_points(d$all_events, grid)
  expect_lt(sum(sm == 0), sum(raw == 0))

  # wider bandwidth spreads support further
  sm_wide <- smooth_offset(d$all_events, grid, bandwidth = 900)
  expect_lte(sum(sm_wide == 0), sum(sm == 0))

  # errors when nothing is within reach
  far_away <- sf::st_sfc(sf::st_point(c(1e6, 1e6)), crs = 32610)
  expect_error(smooth_offset(far_away, grid, bandwidth = 450), "no denominator events")
})

test_that("a smoothed offset can be used to fit an rtm rate model", {
  d <- synthetic_offset_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  sm <- smooth_offset(d$all_events, grid, bandwidth = 450)

  set.seed(5)
  fit <- suppressWarnings(rtm(
    outcome = d$fatal_events,
    offset = sm,
    factors = list(naloxone = list(data = d$naloxone, type = "protective")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  ))

  expect_s3_class(fit, "rtm")
  expect_true(fit$has_offset)
  expect_true("naloxone" %in% fit$selected$factor)
})

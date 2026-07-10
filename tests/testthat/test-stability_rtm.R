fit_synthetic <- function() {
  d <- synthetic_rtm_data(seed = 42)
  set.seed(1)
  rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars, parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
}

test_that("stability_rtm reports frequencies and bagged RRVs", {
  fit <- fit_synthetic()
  set.seed(2)
  st <- stability_rtm(fit, nboot = 10, verbose = FALSE)

  expect_s3_class(st, "rtm_stability")
  expect_setequal(st$factors$factor, c("bars", "parks"))
  expect_true(all(st$factors$freq >= 0 & st$factors$freq <= 1))
  # bars is a strong true effect: selected in (nearly) every resample
  bars <- st$factors[st$factors$factor == "bars", ]
  expect_gte(bars$freq, 0.8)
  expect_gt(bars$rrv_bagged, 1)
  # full-data RRV is carried over for comparison
  expect_equal(bars$rrv_full, fit$selected$rrv[fit$selected$factor == "bars"])
  # percentile interval brackets the bagged estimate
  expect_true(all(st$factors$rrv_lo <= st$factors$rrv_bagged + 1e-12))
  expect_true(all(st$factors$rrv_hi >= st$factors$rrv_bagged - 1e-12))
  # variable-level frequencies aggregate consistently: a factor is selected
  # at most once per resample, so variable freqs sum to at most factor freq
  for (f in st$factors$factor) {
    vsum <- sum(st$variables$freq[st$variables$factor == f])
    expect_lte(vsum, st$factors$freq[st$factors$factor == f] + 1e-12)
  }
  expect_output(print(st), "freq")
})

test_that("stability_rtm supports block resampling", {
  fit <- fit_synthetic()
  set.seed(3)
  st <- stability_rtm(fit, nboot = 5, resample = "block", cluster_size = 600, verbose = FALSE)
  expect_equal(st$resample, "block")
  expect_equal(st$cluster_size, 600)
  expect_true(all(st$factors$freq >= 0 & st$factors$freq <= 1))
})

test_that("stability_rtm works for offset models", {
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
  set.seed(6)
  st <- stability_rtm(fit, nboot = 5, verbose = FALSE)
  nal <- st$factors[st$factors$factor == "naloxone", ]
  expect_gte(nal$freq, 0.6)
  expect_true(is.finite(nal$rrv_bagged))
})

test_that("stability_rtm validates its inputs", {
  fit <- fit_synthetic()
  expect_error(stability_rtm(list()), "rtm object")
  expect_error(stability_rtm(fit, nboot = 1), "at least 2")
  expect_error(stability_rtm(fit, level = 2), "level")
  expect_error(
    stability_rtm(fit, nboot = 2, resample = "block", cluster_size = 1e6, verbose = FALSE),
    "fewer than 2 tiles"
  )
})

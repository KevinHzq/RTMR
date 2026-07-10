fit_pair <- function() {
  d <- synthetic_rtm_data(seed = 42)
  set.seed(1)
  full <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars, parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
  set.seed(1)
  weak <- rtm(
    outcome = d$crimes,
    factors = list(parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
  list(full = full, weak = weak)
}

test_that("check_fit reports pseudo-R2 and calibrated predictive checks", {
  fits <- fit_pair()
  set.seed(2)
  cf <- check_fit(fits$full, nsim = 199)

  expect_s3_class(cf, "rtm_fit_check")
  expect_named(cf$pseudo_r2, c("deviance", "mcfadden", "nagelkerke"))
  expect_true(all(cf$pseudo_r2 > 0 & cf$pseudo_r2 < 1))
  expect_setequal(cf$checks$statistic, c("zero_share", "dispersion", "max_count"))
  expect_true(all(cf$checks$p_value > 0 & cf$checks$p_value <= 1))
  expect_true(all(cf$checks$sim_lo <= cf$checks$sim_hi))

  # the model was fit to the true generating process: its own simulations
  # should reproduce the zero share comfortably
  zero <- cf$checks[cf$checks$statistic == "zero_share", ]
  expect_gt(zero$p_value, 0.01)
  expect_output(print(cf), "Pseudo-R-squared")
})

test_that("a stronger model earns a higher deviance R2 on the same data", {
  fits <- fit_pair()
  set.seed(3)
  r2_full <- check_fit(fits$full, nsim = 9)$pseudo_r2["deviance"]
  set.seed(3)
  r2_weak <- check_fit(fits$weak, nsim = 9)$pseudo_r2["deviance"]
  expect_gt(r2_full, r2_weak)
})

test_that("check_fit flags structural excess zeros", {
  set.seed(7)
  counts <- rpois(100, 2)
  counts[sample.int(100, 40)] <- 0 # structural zeros the model knows nothing about
  x <- fake_rtm(counts)

  set.seed(8)
  cf <- check_fit(x, nsim = 199)
  zero <- cf$checks[cf$checks$statistic == "zero_share", ]
  expect_lte(zero$p_value, 0.05)
  expect_gt(zero$observed, zero$sim_hi)
  expect_output(print(cf), "structural zeros")
})

test_that("check_fit works for offset models on fitted cells only", {
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
  cf <- check_fit(fit, nsim = 99)
  expect_equal(cf$n_cells, sum(fit$grid$offset_count > 0))
  expect_true(is.finite(cf$pseudo_r2["deviance"]))
})

test_that("check_fit validates its inputs", {
  expect_error(check_fit(list()), "rtm object")
  set.seed(9)
  x <- fake_rtm(rpois(100, 2))
  expect_error(check_fit(x, nsim = 1), "nsim")
  expect_error(check_fit(x, level = 0), "level")
})

test_that("expected frequencies match the fitted distribution", {
  mu <- c(0.2, 1.5, 3)
  k <- 0:50
  # Poisson: frequencies over all k sum to the number of cells, and k = 0
  # matches the direct computation
  ef <- RTMR:::expected_frequencies(mu, NULL, k)
  expect_equal(sum(ef), length(mu), tolerance = 1e-8)
  expect_equal(ef[1], sum(dpois(0, mu)))
  # negative binomial ditto
  ef_nb <- RTMR:::expected_frequencies(mu, 1.2, k)
  expect_equal(sum(ef_nb), length(mu), tolerance = 1e-6)
  expect_equal(ef_nb[1], sum(dnbinom(0, mu = mu, size = 1.2)))
})

test_that("quantile residuals are standard normal under a correct model", {
  set.seed(10)
  mu <- rep(2, 2000)
  y <- rpois(2000, mu)
  r <- RTMR:::quantile_residuals(y, mu)
  expect_gt(stats::ks.test(r, "pnorm")$p.value, 0.01)
  expect_equal(mean(r), 0, tolerance = 0.1)
  expect_equal(sd(r), 1, tolerance = 0.1)
})

test_that("plot.rtm_fit_check draws rootogram and QQ plots", {
  set.seed(9)
  x <- fake_rtm(rpois(100, 2))
  set.seed(11)
  cf <- check_fit(x, nsim = 19)

  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(cf))
  expect_silent(plot(cf, which = "rootogram"))
  expect_silent(plot(cf, which = "qq"))
  expect_silent(plot(cf, which = "rootogram", max_count = 10))
  expect_error(plot(cf, which = "nope"), "arg")
})

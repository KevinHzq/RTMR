test_that("moran_rtm matches the brute-force dense-matrix Moran's I", {
  set.seed(3)
  x <- fake_rtm(stats::rpois(100, 2))
  m <- moran_rtm(x, nsim = 9)

  # dense row-standardized queen contiguity weight matrix
  geom <- sf::st_geometry(x$grid)
  nb <- sf::st_intersects(geom)
  n <- length(geom)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, setdiff(nb[[i]], i)] <- 1
  W <- W / rowSums(W)

  e <- stats::residuals(x$best_model, type = "pearson")
  z <- e - mean(e)
  I_dense <- (n / sum(W)) * as.numeric(t(z) %*% W %*% z) / sum(z^2)

  expect_s3_class(m, "rtm_moran")
  expect_equal(m$statistic, I_dense, tolerance = 1e-12)
  expect_equal(m$expected, -1 / (n - 1))
  expect_equal(m$n_cells, n)
  expect_equal(m$n_dropped, 0)
})

test_that("moran_rtm detects strong residual spatial autocorrelation", {
  # counts follow a smooth spatial gradient that the intercept-only model
  # cannot explain, so residuals are strongly clustered
  xy <- expand.grid(col = 1:10, row = 1:10)
  lambda <- exp(0.4 * (xy$col + xy$row) / 2 - 2)
  set.seed(11)
  x <- fake_rtm(stats::rpois(100, lambda))

  set.seed(12)
  m <- moran_rtm(x, nsim = 199)
  expect_gt(m$statistic, m$expected)
  expect_lte(m$p_value, 0.05)
  expect_output(print(m), "Moran's I")
})

test_that("moran_rtm is near null for exchangeable residuals", {
  set.seed(21)
  x <- fake_rtm(stats::rpois(100, 2))
  set.seed(22)
  m <- moran_rtm(x, nsim = 499)
  # i.i.d. counts: no evidence of clustering at a strict level
  expect_gt(m$p_value, 0.01)
})

test_that("moran_rtm works end-to-end on a fitted rtm object", {
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
  set.seed(2)
  m <- moran_rtm(fit, nsim = 99)
  expect_s3_class(m, "rtm_moran")
  expect_true(is.finite(m$statistic))
  expect_true(m$p_value > 0 && m$p_value <= 1)
  expect_equal(m$n_cells + m$n_dropped, nrow(fit$grid))
  expect_length(m$sims, 99)
})

test_that("moran_rtm uses only fitted cells for offset models and drops islands", {
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
  m <- moran_rtm(fit, nsim = 99)
  n_fitted <- sum(fit$grid$offset_count > 0)
  expect_equal(m$n_cells + m$n_dropped, n_fitted)
  expect_lt(m$n_cells, nrow(fit$grid))
})

test_that("moran_rtm validates its inputs", {
  expect_error(moran_rtm(list()), "rtm object")
  set.seed(31)
  x <- fake_rtm(stats::rpois(100, 2))
  expect_error(moran_rtm(x, nsim = 0), "nsim")
  x$grid <- x$grid[1:50, ]
  expect_error(moran_rtm(x), "do not match")
})

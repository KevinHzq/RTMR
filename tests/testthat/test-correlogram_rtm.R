test_that("band-wise Moran's I matches a brute-force dense computation", {
  set.seed(3)
  x <- fake_rtm(stats::rpois(100, 2))
  cg <- correlogram_rtm(x, breaks = c(150, 300), nsim = 9)

  e <- stats::residuals(x$best_model, type = "pearson")
  z <- e - mean(e)
  n <- length(z)
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x$grid))))
  D <- as.matrix(stats::dist(ctr))

  for (k in 1:2) {
    lo <- c(0, 150)[k]
    hi <- c(150, 300)[k]
    W <- (D > lo & D <= hi) * 1
    diag(W) <- 0
    deg <- rowSums(W)
    Wr <- W / ifelse(deg == 0, 1, deg)
    s0 <- sum(Wr)
    I_dense <- (n / s0) * as.numeric(t(z) %*% Wr %*% z) / sum(z^2)
    expect_equal(cg$bands$moran[k], I_dense, tolerance = 1e-12)
    expect_equal(cg$bands$n_pairs[k], sum(W) / 2)
  }
})

test_that("correlogram recovers the range of a short-range residual pattern", {
  # residual correlation confined to 2x2 blocks of cells (200 m): counts
  # share a common block effect, so correlation should die out beyond
  # roughly one block
  xy <- expand.grid(col = 0:9, row = 0:9)
  block <- paste(xy$col %/% 2, xy$row %/% 2)
  set.seed(13)
  block_effect <- stats::rnorm(length(unique(block)), 0, 0.8)[factor(block)]
  counts <- stats::rpois(100, exp(1 + block_effect))
  x <- fake_rtm(counts)

  set.seed(14)
  cg <- correlogram_rtm(x, breaks = c(150, 300, 450, 600, 750), nsim = 199)
  # adjacent cells (first band) are correlated...
  expect_lte(cg$bands$p_value[1], 0.05)
  # ...and the estimated range is finite and modest
  expect_true(is.finite(cg$range) && cg$range > 0)
  expect_lte(cg$range, 450)
  expect_equal(cg$suggested_cluster_size, cg$range)
})

test_that("correlogram reports no range for exchangeable residuals", {
  set.seed(21)
  x <- fake_rtm(stats::rpois(100, 2))
  set.seed(22)
  cg <- correlogram_rtm(x, breaks = c(150, 300, 450), nsim = 199)
  expect_equal(cg$range, 0)
  expect_true(is.na(cg$suggested_cluster_size))
  expect_output(print(cg), "model-based standard errors are likely adequate")
})

test_that("correlogram works end-to-end, prints, and plots", {
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
  cg <- correlogram_rtm(fit, nsim = 49)
  expect_s3_class(cg, "rtm_correlogram")
  # default bands: one cell width apart, up to 3 * max influence (2700 m)
  expect_equal(max(cg$bands$upper), 3 * max(fit$meta$spatial_influence, na.rm = TRUE))
  expect_true(all(diff(cg$bands$upper) > 0))
  expect_true(all(cg$bands$n_pairs > 0))
  expect_output(print(cg), "correlogram")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(cg))
})

test_that("correlogram uses fitted cells only for offset models", {
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
  cg <- correlogram_rtm(fit, breaks = c(300, 600), nsim = 49)
  expect_equal(cg$n_cells, sum(fit$grid$offset_count > 0))
})

test_that("correlogram validates its inputs", {
  expect_error(correlogram_rtm(list()), "rtm object")
  set.seed(31)
  x <- fake_rtm(stats::rpois(100, 2))
  expect_error(correlogram_rtm(x, nsim = 0), "nsim")
  expect_error(correlogram_rtm(x, breaks = c(-100, 200)), "positive")
  expect_error(correlogram_rtm(x, alpha_level = 2), "alpha_level")
})

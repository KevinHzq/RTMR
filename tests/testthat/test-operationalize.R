test_that("spatial_influences builds whole- and half-block sequences", {
  expect_equal(spatial_influences(500, max_blocks = 3), c(500, 1000, 1500))
  expect_equal(
    spatial_influences(500, max_blocks = 3, increment = "half"),
    c(250, 500, 750, 1000, 1250, 1500)
  )
  expect_equal(spatial_influences(100, max_blocks = 4), c(100, 200, 300, 400))
})

test_that("binarize_density flags cells 2 sd above the mean", {
  x <- c(rep(0, 50), rep(1, 49), 100)
  expect_equal(which(binarize_density(x) == 1), 100L)
})

test_that("operationalize builds labelled binary variables with metadata", {
  d <- synthetic_rtm_data()
  grid <- create_grid(d$boundary, cellsize = 300)

  out <- operationalize(
    d$bars, grid, "bars",
    block_length = 300, max_blocks = 2, operation = "both"
  )

  expect_equal(nrow(out), nrow(grid))
  expect_named(
    out,
    c("bars_prox_300", "bars_prox_600", "bars_dens_300", "bars_dens_600")
  )
  expect_true(all(unlist(out) %in% c(0L, 1L)))

  meta <- attr(out, "rtm_meta")
  expect_equal(meta$variable, names(out))
  expect_equal(unique(meta$factor), "bars")
  expect_equal(meta$spatial_influence, c(300, 600, 300, 600))

  # wider proximity thresholds cover at least as many cells
  expect_true(sum(out$bars_prox_600) >= sum(out$bars_prox_300))
})

test_that("polygon factors work with proximity and are rejected for density", {
  d <- synthetic_rtm_data()
  grid <- create_grid(d$boundary, cellsize = 300)

  # two square polygon features
  sq <- function(x0, y0, w) {
    sf::st_polygon(list(rbind(
      c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0)
    )))
  }
  polys <- sf::st_sfc(sq(300, 300, 600), sq(1800, 1800, 450), crs = 32610)

  out <- operationalize(
    polys, grid, "parks",
    block_length = 300, max_blocks = 2, operation = "proximity"
  )
  expect_true(all(unlist(out) %in% c(0L, 1L)))

  # cells whose centroid falls inside a polygon are at distance 0: exposed
  # at every threshold
  centroids <- suppressWarnings(sf::st_centroid(sf::st_geometry(grid)))
  inside <- lengths(sf::st_intersects(centroids, polys)) > 0
  expect_true(any(inside))
  expect_true(all(out$parks_prox_300[inside] == 1L))
  expect_true(all(out$parks_prox_600[inside] == 1L))

  # density (and both) are refused for non-point geometry
  expect_error(
    operationalize(polys, grid, "parks", block_length = 300, operation = "density"),
    "requires point features"
  )
  expect_error(
    operationalize(polys, grid, "parks", block_length = 300, operation = "both"),
    "requires point features"
  )
})

test_that("protective factors are coded 0/-1", {
  d <- synthetic_rtm_data()
  grid <- create_grid(d$boundary, cellsize = 300)

  out <- operationalize(
    d$parks, grid, "parks",
    block_length = 300, max_blocks = 1, type = "protective"
  )
  expect_true(all(unlist(out) %in% c(0L, -1L)))
  expect_true(any(unlist(out) == -1L))
})

test_that("compute_density matches a direct Epanechnikov kernel sum", {
  d <- synthetic_rtm_data()
  grid <- create_grid(d$boundary, cellsize = 300)
  h <- 500

  dens <- compute_density(d$bars, grid, bandwidth = h)

  centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(grid))))
  pts <- sf::st_coordinates(d$bars)
  expected <- apply(centroids, 1, function(cc) {
    u <- sqrt((pts[, 1] - cc[1])^2 + (pts[, 2] - cc[2])^2) / h
    sum(2 / (pi * h^2) * (1 - u[u <= 1]^2))
  })
  expect_equal(dens, expected, tolerance = 1e-8)
})

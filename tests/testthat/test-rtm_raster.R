test_that("rtm_raster reproduces the grid values as a raster", {
  skip_if_not_installed("terra")

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

  r <- rtm_raster(fit)
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::res(r), c(150, 150))

  # every grid cell value must appear at its centroid in the raster
  centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(fit$grid))))
  vals <- terra::extract(r, centroids)[, 1]
  expect_equal(vals, fit$grid$relrisk, tolerance = 1e-6)

  # prediction layer too
  r_pred <- rtm_raster(fit, what = "prediction")
  vals_pred <- terra::extract(r_pred, centroids)[, 1]
  expect_equal(vals_pred, fit$grid$prediction, tolerance = 1e-6)
})

test_that("write_rtm writes a readable GeoTiff", {
  skip_if_not_installed("terra")

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

  tif <- file.path(tempdir(), "rtm_test_score.tif")
  on.exit(unlink(tif))

  write_rtm(fit, tif, overwrite = TRUE)
  expect_true(file.exists(tif))

  r <- terra::rast(tif)
  expect_equal(terra::res(r), c(150, 150))
  expect_equal(
    max(terra::values(r), na.rm = TRUE),
    max(fit$grid$relrisk),
    tolerance = 1e-6
  )
})

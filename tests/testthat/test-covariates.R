# crimes driven by bars AND a north-south gradient (confounder surface),
# plus extra events inside a binary "zone" (candidate covariate)
synthetic_covariate_data <- function(seed = 11) {
  set.seed(seed)
  size <- 3000

  boundary <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0, 0), c(size, 0), c(size, size), c(0, size), c(0, 0)
    ))),
    crs = 32610
  ))

  rand_pts <- function(n) cbind(stats::runif(n, 0, size), stats::runif(n, 0, size))
  as_sfc <- function(m) {
    sf::st_sfc(lapply(seq_len(nrow(m)), function(i) sf::st_point(m[i, ])), crs = 32610)
  }

  bars <- rand_pts(15)

  crimes <- do.call(rbind, lapply(seq_len(nrow(bars)), function(i) {
    n <- stats::rpois(1, 8)
    cbind(stats::rnorm(n, bars[i, 1], 100), stats::rnorm(n, bars[i, 2], 100))
  }))
  # north-south gradient: background events concentrated in the north
  north_bg <- cbind(stats::runif(80, 0, size), stats::runif(80, size / 2, size))
  # extra events inside the south-west quadrant "zone"
  zone_ev <- cbind(stats::runif(40, 0, size / 2), stats::runif(40, 0, size / 2))
  crimes <- rbind(crimes, north_bg, zone_ev)
  crimes <- crimes[
    crimes[, 1] >= 0 & crimes[, 1] <= size & crimes[, 2] >= 0 & crimes[, 2] <= size, ,
    drop = FALSE
  ]

  list(boundary = boundary, bars = as_sfc(bars), crimes = as_sfc(crimes))
}

test_that("adjustment covariates are always included but never selected", {
  d <- synthetic_covariate_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(grid))))
  northing <- as.numeric(scale(centroids[, 2]))

  set.seed(2)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars),
    covariates = list(northing = list(values = northing, role = "adjustment")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )

  # reported in adjustments with the planted positive direction, not in the
  # RRV table
  expect_equal(fit$adjustments$covariate, "northing")
  expect_gt(fit$adjustments$coefficient, 0)
  expect_false("northing" %in% fit$selected$factor)
  expect_true("bars" %in% fit$selected$factor)

  # present in the fitted model itself
  expect_true("northing" %in% names(coef(fit$best_model)))

  expect_output(print(fit), "Adjustment covariates")
  expect_output(print(fit), "northing")
})

test_that("binary candidate covariates compete for selection and get an RRV", {
  d <- synthetic_covariate_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(grid))))
  zone <- as.integer(centroids[, 1] < 1500 & centroids[, 2] < 1500)

  set.seed(2)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars),
    covariates = list(zone = list(values = zone, role = "candidate")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )

  expect_true("zone" %in% fit$selected$factor)
  zrow <- fit$selected[fit$selected$factor == "zone", ]
  expect_equal(zrow$operation, "covariate")
  expect_true(is.na(zrow$spatial_influence))
  expect_gt(zrow$rrv, 1)
  expect_null(fit$adjustments)
})

test_that("covariate input validation and warnings", {
  d <- synthetic_covariate_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  n <- nrow(grid)

  base_args <- list(
    outcome = d$crimes, factors = list(bars = d$bars),
    boundary = d$boundary, cell_size = 150, block_length = 300,
    operation = "proximity", verbose = FALSE
  )

  # bare vector: role must be specified explicitly
  expect_error(
    do.call(rtm, c(base_args, list(covariates = list(z = rep(1, n))))),
    "`values` and `role`"
  )
  # wrong length
  expect_error(
    do.call(rtm, c(base_args, list(
      covariates = list(z = list(values = 1:3, role = "adjustment"))
    ))),
    "one value per grid cell"
  )
  # continuous candidate warns about per-unit RRV
  expect_warning(
    do.call(rtm, c(base_args, list(
      covariates = list(z = list(values = runif(n), role = "candidate"))
    ))),
    "not binary"
  )
  # constant adjustment is dropped with a warning
  expect_warning(
    do.call(rtm, c(base_args, list(
      covariates = list(z = list(values = rep(1, n), role = "adjustment"))
    ))),
    "constant adjustment"
  )
})

test_that("categorical adjustments are dummy-coded against the chosen reference", {
  d <- synthetic_covariate_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  xy <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(grid))))

  # three land-use classes; the north (where background crime was planted)
  # is commercial
  landuse <- ifelse(xy[, 2] > 1500, "commercial",
    ifelse(xy[, 1] > 1500, "industrial", "residential")
  )

  base_args <- list(
    outcome = d$crimes, factors = list(bars = d$bars),
    boundary = d$boundary, cell_size = 150, block_length = 300,
    operation = "proximity", verbose = FALSE
  )

  set.seed(2)
  fit <- do.call(rtm, c(base_args, list(covariates = list(
    landuse = list(values = landuse, role = "adjustment", reference = "residential")
  ))))

  # one dummy per non-reference level, none for the reference itself
  expect_setequal(
    fit$adjustments$covariate,
    c("landuse_commercial", "landuse_industrial")
  )
  # industrial (no background events planted there) has a negative log rate
  # ratio vs residential (which received the zone events)
  expect_lt(
    fit$adjustments$coefficient[fit$adjustments$covariate == "landuse_industrial"],
    0
  )
  # dummies are in the fitted model, bars still selected
  expect_true(all(c("landuse_commercial", "landuse_industrial") %in%
    names(coef(fit$best_model))))
  expect_true("bars" %in% fit$selected$factor)

  # default reference is the first level (alphabetical: commercial)
  set.seed(2)
  fit_default <- do.call(rtm, c(base_args, list(covariates = list(
    landuse = list(values = landuse, role = "adjustment")
  ))))
  expect_setequal(
    fit_default$adjustments$covariate,
    c("landuse_industrial", "landuse_residential")
  )

  # unknown reference level errors, listing the available levels
  expect_error(
    do.call(rtm, c(base_args, list(covariates = list(
      landuse = list(values = landuse, role = "adjustment", reference = "parkland")
    )))),
    "not found.*commercial, industrial, residential"
  )

  # categorical candidates are refused with a pointer to binary indicators
  expect_error(
    do.call(rtm, c(base_args, list(covariates = list(
      landuse = list(values = landuse, role = "candidate")
    )))),
    "binary indicator"
  )

  # missing values are refused
  landuse_na <- landuse
  landuse_na[1] <- NA
  expect_error(
    do.call(rtm, c(base_args, list(covariates = list(
      landuse = list(values = landuse_na, role = "adjustment")
    )))),
    "missing values"
  )

  # reference on a numeric covariate warns that it is ignored
  expect_warning(
    do.call(rtm, c(base_args, list(covariates = list(
      north = list(values = as.numeric(scale(xy[, 2])), role = "adjustment", reference = "residential")
    )))),
    "ignored"
  )
})

test_that("protective candidate covariates are negated like protective factors", {
  d <- synthetic_covariate_data()
  grid <- create_grid(d$boundary, cellsize = 150)
  zone <- as.integer(sf::st_coordinates(
    suppressWarnings(sf::st_centroid(sf::st_geometry(grid)))
  )[, 1] < 1500)

  set.seed(2)
  fit <- rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars),
    covariates = list(zone = list(values = zone, role = "candidate", type = "protective")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )

  # coded 0/-1 in the model data regardless of selection
  expect_true(all(fit$data$zone %in% c(0L, -1L)))
  expect_equal(fit$meta$type[fit$meta$variable == "zone"], "protective")
})

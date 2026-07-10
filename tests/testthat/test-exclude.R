square <- function(x0, y0, x1, y1, crs = 32610) {
  sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)
  ))), crs = crs)
}

test_that("create_grid drops cells fully covered by the exclusion only", {
  boundary <- sf::st_sf(geometry = square(0, 0, 1000, 1000))
  lake <- square(0, 0, 525, 525) # covers 3x3 cells fully, partial beyond

  full <- create_grid(boundary, cellsize = 100)
  masked <- create_grid(boundary, cellsize = 100, exclude = lake)

  expect_equal(nrow(full), 100)
  # exactly the 5x5 block of cells inside (0, 500)^2 is fully covered
  expect_equal(nrow(masked), 100 - 25)
  # no remaining cell is fully covered ...
  expect_true(all(lengths(sf::st_covered_by(masked, lake)) == 0))
  # ... but partially covered cells are kept
  expect_gt(sum(lengths(sf::st_intersects(masked, lake)) > 0), 0)
  # cell ids are preserved from the unmasked grid
  expect_true(all(masked$cell_id %in% full$cell_id))
})

test_that("create_grid errors when the exclusion removes everything", {
  boundary <- sf::st_sf(geometry = square(0, 0, 1000, 1000))
  everything <- square(-100, -100, 1100, 1100)
  expect_error(
    create_grid(boundary, cellsize = 100, exclude = everything),
    "covers every grid cell"
  )
  expect_error(
    create_grid(boundary, cellsize = 100, exclude = 1:3),
    "sf or sfc"
  )
})

test_that("rtm excludes masked cells and warns about events inside the mask", {
  d <- synthetic_rtm_data(seed = 42)
  lake <- square(0, 0, 610, 610)

  # place an event squarely inside the lake to trigger the check
  crimes <- c(d$crimes, sf::st_sfc(sf::st_point(c(100, 100)), crs = 32610))

  expect_warning(
    {
      set.seed(1)
      fit <- rtm(
        outcome = crimes,
        factors = list(bars = d$bars, parks = d$parks),
        boundary = d$boundary,
        cell_size = 150,
        block_length = 300,
        operation = "proximity",
        exclude = lake,
        verbose = FALSE
      )
    },
    "inside the excluded area"
  )

  # the 4x4 block of 150 m cells within (0, 600)^2 is gone
  expect_equal(nrow(fit$grid), 400 - 16)
  expect_true(all(lengths(sf::st_covered_by(fit$grid, lake)) == 0))
  # outputs stay aligned with the reduced grid
  expect_length(fit$grid$relrisk, 400 - 16)
  expect_equal(nrow(fit$data), 400 - 16)

  # downstream diagnostics run on the masked grid
  set.seed(2)
  m <- suppressMessages(moran_rtm(fit, nsim = 49))
  expect_equal(m$n_cells + m$n_dropped, 400 - 16)
})

test_that("rtm with exclude does not warn when no events touch the mask", {
  d <- synthetic_rtm_data(seed = 42)
  # a corner the seed-42 events do not reach
  corner <- square(2850, 0, 3000, 150)
  inside <- lengths(sf::st_intersects(sf::st_geometry(d$crimes), corner)) > 0
  skip_if(any(inside), "seed places an event in the test corner")

  set.seed(1)
  expect_no_warning(
    rtm(
      outcome = d$crimes,
      factors = list(bars = d$bars, parks = d$parks),
      boundary = d$boundary,
      cell_size = 150,
      block_length = 300,
      operation = "proximity",
      exclude = corner,
      verbose = FALSE
    )
  )
})

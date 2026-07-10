test_that("count_points counts each point exactly once", {
  b <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)))),
    crs = 32610
  ))
  grid <- create_grid(b, cellsize = 1)

  pts <- sf::st_sfc(
    sf::st_point(c(1, 0.5)), # exactly on the edge shared by two cells
    sf::st_point(c(1, 1)), # on the corner shared by all four cells
    sf::st_point(c(0.5, 0.5)), # interior
    sf::st_point(c(10, 10)), # outside the grid: dropped
    crs = 32610
  )

  cnt <- count_points(pts, grid)
  expect_length(cnt, nrow(grid))
  expect_equal(sum(cnt), 3)

  # edge assignment is deterministic (first intersecting cell in grid order)
  expect_equal(cnt, count_points(pts, grid))
})

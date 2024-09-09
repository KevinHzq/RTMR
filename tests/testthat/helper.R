# load helper data

test_polygon <- function() {
  sf::st_read(system.file("shape/nc.shp", package="sf"))
}

test_point <- function(x, n = 100, ...) {
  sf::st_sample(x, size = n, ...)
}

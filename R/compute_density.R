# count number of points in each grid
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r

count_points <- function(pt, grid) {
  lengths(sf::st_intersects(grid, pt))
}

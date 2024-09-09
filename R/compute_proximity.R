compute_proximity <- function(pt, grid) {
  grid_centroid <- st_centroid(st_geometry(grid))

  nearest_pt <- sf::st_nearest_feature(grid_centroid, pt)

  sf::st_distance(grid_centroid, pt[nearest_pt], by_element = TRUE)
}

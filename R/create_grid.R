creat_grid <- function(x, ..., clip = TRUE) {
  xboundary <- sf::st_union(x)

  sfc_grid <- sf::st_make_grid(xboundary, ...)

  sf_grid <- sf::st_sf(cell_id = 1:length(sfc_grid), geometry = sfc_grid)

  if (clip) {
    sf::st_filter(sf_grid, xboundary)
  } else
    return(sf_grid)
}

# load helper data

test_polygon <- function() {
  sf::st_read(system.file("shape/nc.shp", package="sf"))
}

test_point <- function(x, n = 100, ...) {
  sf::st_sample(x, size = n, ...)
}

# build a minimal rtm-like object from a grid and per-cell counts, fitted
# with an intercept-only Poisson model, so residual structure is controlled
fake_rtm <- function(counts, nx = 10, ny = 10, cell = 100) {
  boundary <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0, 0), c(nx * cell, 0), c(nx * cell, ny * cell), c(0, ny * cell), c(0, 0)
    ))),
    crs = 32610
  ))
  grid <- create_grid(boundary, cellsize = cell)
  stopifnot(nrow(grid) == length(counts))
  grid$outcome_count <- counts
  fit <- stats::glm(counts ~ 1, family = stats::poisson())
  structure(
    list(best_model = fit, grid = grid, has_offset = FALSE),
    class = "rtm"
  )
}

# synthetic study area (3 km square, metric CRS) with a risk factor ("bars")
# that generates most outcome events nearby, plus an irrelevant factor
# ("parks") and background noise
synthetic_rtm_data <- function(seed = 42,
                               n_bars = 15,
                               n_parks = 10,
                               events_per_bar = 8,
                               sd_around_bar = 100,
                               n_background = 30) {
  set.seed(seed)

  size <- 3000
  boundary <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0, 0), c(size, 0), c(size, size), c(0, size), c(0, 0)
    ))),
    crs = 32610
  ))

  rand_pts <- function(n) cbind(stats::runif(n, 0, size), stats::runif(n, 0, size))
  as_sfc <- function(m) sf::st_sfc(lapply(seq_len(nrow(m)), function(i) sf::st_point(m[i, ])), crs = 32610)

  bars <- rand_pts(n_bars)
  parks <- rand_pts(n_parks)

  n_events <- stats::rpois(n_bars, events_per_bar)
  crimes <- do.call(rbind, lapply(seq_len(n_bars), function(i) {
    cbind(
      stats::rnorm(n_events[i], bars[i, 1], sd_around_bar),
      stats::rnorm(n_events[i], bars[i, 2], sd_around_bar)
    )
  }))
  crimes <- rbind(crimes, rand_pts(n_background))
  crimes <- crimes[
    crimes[, 1] >= 0 & crimes[, 1] <= size & crimes[, 2] >= 0 & crimes[, 2] <= size, ,
    drop = FALSE
  ]

  list(
    boundary = boundary,
    bars = as_sfc(bars),
    parks = as_sfc(parks),
    crimes = as_sfc(crimes)
  )
}

# case-fatality scenario: overdose events cluster around "bars"; an event's
# probability of being fatal is much higher far from a naloxone site, so the
# lethality (rate) model should select naloxone as protective even though
# event volume is driven by bars
synthetic_offset_data <- function(seed = 7) {
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
  naloxone <- rand_pts(8)

  gen_events <- function() {
    events <- do.call(rbind, lapply(seq_len(nrow(bars)), function(i) {
      n <- stats::rpois(1, 40)
      cbind(stats::rnorm(n, bars[i, 1], 150), stats::rnorm(n, bars[i, 2], 150))
    }))
    events <- rbind(events, rand_pts(150))
    events[
      events[, 1] >= 0 & events[, 1] <= size & events[, 2] >= 0 & events[, 2] <= size, ,
      drop = FALSE
    ]
  }

  # fatality probability depends on distance to nearest naloxone site
  split_fatal <- function(events) {
    d_nalox <- apply(events, 1, function(e) {
      min(sqrt((naloxone[, 1] - e[1])^2 + (naloxone[, 2] - e[2])^2))
    })
    p_fatal <- ifelse(d_nalox <= 600, 0.03, 0.35)
    stats::runif(nrow(events)) < p_fatal
  }

  events <- gen_events()
  fatal <- split_fatal(events)

  # independent second period from the same generating process, for holdout
  # validation
  events2 <- gen_events()
  fatal2 <- split_fatal(events2)

  list(
    boundary = boundary,
    bars = as_sfc(bars),
    naloxone = as_sfc(naloxone),
    all_events = as_sfc(events),
    fatal_events = as_sfc(events[fatal, , drop = FALSE]),
    all_events2 = as_sfc(events2),
    fatal_events2 = as_sfc(events2[fatal2, , drop = FALSE])
  )
}

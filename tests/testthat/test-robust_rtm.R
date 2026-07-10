skip_if_not_installed("sandwich")

fit_synthetic <- function() {
  d <- synthetic_rtm_data(seed = 42)
  set.seed(1)
  rtm(
    outcome = d$crimes,
    factors = list(bars = d$bars, parks = d$parks),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
}

test_that("robust_rtm reports robust intervals without changing estimates", {
  fit <- fit_synthetic()
  rob <- suppressWarnings(robust_rtm(fit))

  expect_s3_class(rob, "rtm_robust")
  # point estimates and RRVs are untouched
  expect_equal(rob$selected$rrv, fit$selected$rrv)
  expect_equal(rob$selected$coefficient, fit$selected$coefficient)
  # intervals bracket the RRV and p-values are proper
  expect_true(all(rob$selected$ci_lower <= rob$selected$rrv))
  expect_true(all(rob$selected$ci_upper >= rob$selected$rrv))
  expect_true(all(rob$selected$se_robust > 0))
  expect_true(all(rob$selected$p_robust > 0 & rob$selected$p_robust <= 1))
  # default cluster size is twice the largest candidate influence
  expect_equal(rob$cluster_size, 2 * max(fit$meta$spatial_influence, na.rm = TRUE))
  expect_output(print(rob), "robust")
})

test_that("robust_rtm matches a direct sandwich::vcovCL computation", {
  fit <- fit_synthetic()
  # custom clusters: four quadrants of the study area
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(fit$grid))))
  cl <- paste(ctr[, 1] > 1500, ctr[, 2] > 1500)
  rob <- suppressWarnings(robust_rtm(fit, cluster = cl))

  vc <- sandwich::vcovCL(fit$best_model, cluster = factor(cl))
  se <- sqrt(diag(vc))[rob$selected$variable]
  expect_equal(rob$selected$se_robust, unname(se), tolerance = 1e-12)
  expect_equal(rob$n_clusters, 4)
  expect_true(is.na(rob$cluster_size))
})

test_that("robust_rtm warns on few clusters and errors on degenerate input", {
  fit <- fit_synthetic()
  expect_warning(robust_rtm(fit, cluster_size = 1500), "few clusters")
  expect_error(
    suppressWarnings(robust_rtm(fit, cluster_size = 1e6)),
    "fewer than 2 clusters"
  )
  expect_error(robust_rtm(fit, cluster = 1:5), "one value per grid cell")
  expect_error(robust_rtm(fit, level = 1.5), "level")
  expect_error(robust_rtm(list()), "rtm object")
})

test_that("robust_rtm works for offset models and clusters fitted cells only", {
  d <- synthetic_offset_data(seed = 7)
  set.seed(5)
  fit <- rtm(
    outcome = d$fatal_events,
    offset = d$all_events,
    factors = list(bars = d$bars, naloxone = list(data = d$naloxone, type = "protective")),
    boundary = d$boundary,
    cell_size = 150,
    block_length = 300,
    operation = "proximity",
    verbose = FALSE
  )
  rob <- suppressWarnings(robust_rtm(fit, cluster_size = 600))
  expect_equal(rob$selected$rrv, fit$selected$rrv)
  expect_true(all(is.finite(rob$selected$se_robust)))
  # clusters counted among fitted cells only: cannot exceed tiles overall
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(fit$grid))))
  clusters_all <- length(unique(paste(
    floor((ctr[, 1] - min(ctr[, 1])) / 600),
    floor((ctr[, 2] - min(ctr[, 2])) / 600)
  )))
  expect_lte(rob$n_clusters, clusters_all)
})

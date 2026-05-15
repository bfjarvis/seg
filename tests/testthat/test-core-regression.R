toy_grid <- function(cols = 1:2) {
  data(segdata)
  list(
    coords = as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5))),
    values = as.matrix(segdata[, cols])
  )
}

test_that("localenv accepts coordinate matrices and preserves dimensions", {
  toy <- toy_grid()

  env <- localenv(toy$coords, toy$values, weighting = "unweighted",
                  normalize = FALSE)

  expect_s4_class(env, "SegLocal")
  expect_equal(dim(env@coords), c(100L, 2L))
  expect_equal(dim(env@data), dim(toy$values))
  expect_equal(dim(env@env), dim(toy$values))
  expected <- matrix(colMeans(toy$values), nrow = 100, ncol = 2, byrow = TRUE)
  dimnames(expected) <- dimnames(env@env)
  expect_equal(env@env, expected)
})

test_that("automatic bandwidth is deterministic and drives default localenv", {
  toy <- toy_grid()
  bandwidth <- exp(mean(log(as.numeric(dist(toy$coords)))))
  coords <- cbind(seq_len(150), sin(seq_len(150)))

  env_default <- localenv(toy$coords, toy$values)
  env_explicit <- localenv(toy$coords, toy$values, maxdist = bandwidth)
  set.seed(42)
  old_seed <- .Random.seed
  bandwidth_1 <- seg:::.geometric_mean_distance(coords)
  seed_after <- .Random.seed
  bandwidth_2 <- seg:::.geometric_mean_distance(coords)

  expect_equal(env_default@env, env_explicit@env, tolerance = 1e-12)
  expect_equal(seed_after, old_seed)
  expect_equal(bandwidth_2, bandwidth_1)
  expect_true(is.finite(bandwidth_1))
  expect_gt(bandwidth_1, 0)
})

test_that("default localenv is invariant to coordinate units", {
  toy <- toy_grid()

  env_m <- localenv(toy$coords, toy$values)
  env_km <- localenv(toy$coords * 1000, toy$values)

  expect_equal(env_km@env, env_m@env, tolerance = 1e-12)
})

test_that("maxdist zero uses observed data as local environments", {
  toy <- toy_grid()

  env <- localenv(toy$coords, toy$values, maxdist = 0)
  seg <- spseg(toy$coords, toy$values, maxdist = 0)

  expect_equal(env@env, toy$values)
  expect_equal(seg@env, toy$values)
})

test_that("localenv weighting schemes match between coordinate and dist paths", {
  coords <- cbind(x = c(0, 1, 3), y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     10, 0),
                   ncol = 2, byrow = TRUE)
  colnames(values) <- c("a", "b")

  unweighted <- localenv(coords, values, weighting = "unweighted", maxdist = 1)
  expect_equal(unweighted@env[1, ], c(a = 5, b = 5))
  expect_equal(unweighted@env[3, ], c(a = 10, b = 0))

  for (scheme in c("biweight", "inverse", "exponential", "unweighted")) {
    env_coords <- localenv(coords, values, weighting = scheme, maxdist = 2,
                           power = 2)
    env_dist <- localenv(coords, values, sprel = dist(coords),
                         weighting = scheme, maxdist = 2, power = 2)
    expect_equal(env_dist@env, env_coords@env, tolerance = 1e-12)
  }
})

test_that("localenv ignores deprecated weighting syntax with messages", {
  coords <- cbind(x = c(0, 1, 3), y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     10, 0),
                   ncol = 2, byrow = TRUE)

  expect_message(
    old_biweight <- localenv(coords, values, useExp = FALSE, scale = TRUE,
                             maxdist = 2),
    "Ignoring deprecated localenv argument"
  )
  expect_message(
    old_inverse <- localenv(coords, values, useExp = FALSE, scale = FALSE,
                            maxdist = 2),
    "Ignoring deprecated localenv argument"
  )

  default <- localenv(coords, values, maxdist = 2)
  expect_equal(old_biweight@env,
               default@env)
  expect_equal(old_inverse@env,
               default@env)
})

test_that("localenv and spseg vary across bandwidths", {
  toy <- toy_grid()

  env_small <- localenv(toy$coords, toy$values, maxdist = 1)
  env_large <- localenv(toy$coords, toy$values, maxdist = 5)
  seg_small <- as.list(spseg(env_small))
  seg_large <- as.list(spseg(env_large))

  expect_false(isTRUE(all.equal(env_small@env, env_large@env)))
  expect_false(isTRUE(all.equal(seg_small$d, seg_large$d)))
})

test_that("spseg selected measures match all-measure C results", {
  toy <- toy_grid()
  env <- localenv(toy$coords, toy$values, maxdist = 2)

  all_measures <- as.list(spseg(env, useC = TRUE))
  expect_equal(as.list(spseg(env, measures = "exposure", useC = TRUE))$p,
               all_measures$p, tolerance = 1e-12)
  expect_equal(as.list(spseg(env, measures = "information", useC = TRUE))$h,
               all_measures$h, tolerance = 1e-12)
  expect_equal(as.list(spseg(env, measures = "diversity", useC = TRUE))$r,
               all_measures$r, tolerance = 1e-12)
  expect_equal(as.list(spseg(env, measures = "dissimilarity", useC = TRUE))$d,
               all_measures$d, tolerance = 1e-12)
})

test_that("spseg handles zero entropy terms without tolerance adjustment", {
  coords <- cbind(x = 1:4, y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     5, 0,
                     0, 5),
                   ncol = 2, byrow = TRUE)
  colnames(values) <- c("a", "b")
  env <- SegLocal(coords, values, values)

  c_result <- as.list(spseg(env, measures = "information", useC = TRUE))
  r_result <- as.list(spseg(env, measures = "information", useC = FALSE))

  expect_true(is.finite(c_result$h))
  expect_equal(c_result$h, r_result$h, tolerance = 1e-12)
  expect_equal(c_result$h, 1, tolerance = 1e-12)
})

test_that("spseg dissimilarity uses generalized local-environment formula", {
  coords <- cbind(x = 1:3, y = 0)
  values <- matrix(c(10, 20, 70,
                     15, 15, 70,
                     30, 10, 60),
                   ncol = 3, byrow = TRUE)
  colnames(values) <- c("a", "b", "c")
  env_values <- matrix(c(20, 10, 70,
                         15, 15, 70,
                         10, 20, 70),
                       ncol = 3, byrow = TRUE)
  colnames(env_values) <- colnames(values)
  env <- SegLocal(coords, values, env_values)
  env_prop <- env_values / rowSums(env_values)
  overall_prop <- colSums(values) / sum(values)
  interaction <- sum(overall_prop * (1 - overall_prop))
  expected <- sum(colSums(
    abs(env_prop - matrix(overall_prop, nrow = nrow(env_prop),
                          ncol = ncol(env_prop), byrow = TRUE)) *
      (rowSums(values) / (2 * sum(values) * interaction))
  ))

  expect_equal(as.list(spseg(env, measures = "dissimilarity",
                             useC = FALSE))$d,
               expected, tolerance = 1e-12)
  expect_equal(as.list(spseg(env, measures = "dissimilarity",
                             useC = TRUE))$d,
               expected, tolerance = 1e-12)
  expect_equal(
    spseg(coords, data = values, bands = 1, output = "indices",
          measures = "dissimilarity")@indices$overall$d[[1]],
    as.list(spseg(localenv(coords, data = values, maxdist = 1),
                  measures = "dissimilarity"))$d,
    tolerance = 1e-12
  )
})

test_that("pairwise scalar indices agree with overall indices for two groups", {
  toy <- toy_grid(cols = 5:6)
  result <- spseg(toy$coords, toy$values, bands = c(1, 2),
                  output = "indices", comparison = "both")

  for (band in names(result@indices$overall$d)) {
    expect_equal(result@indices$pairwise$d[[band]][1, 2],
                 result@indices$overall$d[[band]], tolerance = 1e-12)
    expect_equal(result@indices$pairwise$r[[band]][1, 2],
                 result@indices$overall$r[[band]], tolerance = 1e-12)
    expect_equal(result@indices$pairwise$h[[band]][1, 2],
                 result@indices$overall$h[[band]], tolerance = 1e-12)
  }
})

test_that("pairwise indices skip units outside the pair and mark empty pairs", {
  coords <- cbind(x = 1:4, y = 0)
  values <- matrix(c(0, 0, 10, 0,
                     5, 0, 0, 0,
                     0, 5, 0, 0,
                     0, 0, 5, 0),
                   ncol = 4, byrow = TRUE)
  colnames(values) <- c("a", "b", "c", "empty")

  result <- spseg(coords, values, bands = 0, output = "indices",
                  comparison = "pairwise",
                  measures = c("dissimilarity", "information", "diversity"))

  expect_equal(result@indices$pairwise$d[["0"]]["a", "b"], 1,
               tolerance = 1e-12)
  expect_true(is.na(result@indices$pairwise$d[["0"]]["a", "empty"]))
  expect_true(is.finite(result@indices$pairwise$r[["0"]]["a", "b"]))
  expect_true(is.finite(result@indices$pairwise$h[["0"]]["a", "b"]))
  expect_equal(nrow(as.data.frame(result, what = "pairwise")), 18L)
})

test_that("spseg ignores deprecated tolerance argument", {
  coords <- cbind(x = 1:4, y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     5, 0,
                     0, 5),
                   ncol = 2, byrow = TRUE)
  env <- SegLocal(coords, values, values)

  expect_message(
    with_tol <- spseg(env, measures = "information", tol = 1e-3),
    "ignored"
  )
  without_tol <- spseg(env, measures = "information")

  expect_equal(as.list(with_tol)$h, as.list(without_tol)$h)
})

test_that("spseg accommodates deprecated method argument", {
  toy <- toy_grid()
  env <- localenv(toy$coords, toy$values, maxdist = 2)

  expect_warning(
    old <- spseg(env, method = "information", useC = TRUE),
    "deprecated"
  )
  new <- spseg(env, measures = "information", useC = TRUE)

  expect_equal(as.list(old)$h, as.list(new)$h, tolerance = 1e-12)
})

test_that("spatseg remains as a deprecated compatibility wrapper", {
  toy <- toy_grid()
  env <- localenv(toy$coords, toy$values, maxdist = 2)

  expect_warning(
    result <- spatseg(env, method = "dissimilarity", useC = TRUE),
    "deprecated"
  )
  expect_equal(as.list(result)$d,
               as.list(spseg(env, measures = "dissimilarity", useC = TRUE))$d,
               tolerance = 1e-12)
})

test_that("spseg wrapper returns selected measures only", {
  toy <- toy_grid()

  result <- spseg(toy$coords, toy$values,
                  measures = c("information", "dissimilarity"), maxdist = 2)

  expect_s4_class(result, "SegSpatial")
  expect_length(result@h, 1L)
  expect_length(result@d, 1L)
  expect_length(result@r, 0L)
  expect_equal(dim(result@p), c(0L, 0L))
})

test_that("spseg argument order keeps outputs before local environment settings", {
  expect_equal(
    names(formals(spseg))[1:9],
    c("x", "data", "measures", "sprel", "maxdist", "bands", "weighting",
      "power", "normalize")
  )
})

test_that("spseg exposes and syncs localenv weighting arguments", {
  toy <- toy_grid()
  env <- localenv(toy$coords, toy$values, maxdist = 2, weighting = "inverse",
                  power = 2, normalize = FALSE)
  wrapped <- spseg(toy$coords, toy$values, maxdist = 2,
                   weighting = "inverse", power = 2, normalize = FALSE)

  expect_equal(wrapped@env, env@env, tolerance = 1e-12)
})

test_that("spseg unified result can store local environments and indices", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values, bands = c(1, 2),
                  output = "full")
  env2 <- localenv(toy$coords, toy$values, maxdist = 2)
  legacy2 <- as.list(spseg(env2))

  expect_s4_class(result, "SegResult")
  expect_equal(result@bands, c(1, 2))
  expect_length(result@env, 2L)
  expect_equal(result@env[[2]], env2@env, tolerance = 1e-12)
  expect_equal(result@indices$overall$d[["2"]], legacy2$d, tolerance = 1e-12)
  expect_equal(result@indices$overall$r[["2"]], legacy2$r, tolerance = 1e-12)
  expect_equal(result@indices$overall$h[["2"]], legacy2$h, tolerance = 1e-12)
  expect_equal(result@indices$overall$p[["2"]], legacy2$p, tolerance = 1e-12)
})

test_that("spseg indices output avoids storing inputs and environments", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values, bands = c(1, 2),
                  output = "indices",
                  measures = c("information", "dissimilarity"))

  expect_s4_class(result, "SegResult")
  expect_null(result@coords)
  expect_null(result@data)
  expect_null(result@env)
  expect_length(result@indices$overall$h, 2L)
  expect_length(result@indices$overall$d, 2L)
  expect_length(result@indices$overall$r, 0L)
  expect_length(result@indices$overall$p, 0L)
  expect_equal(
    as.data.frame(result)$measure,
    rep(c("d", "h"), each = 2)
  )
})

test_that("default bands use sampled distance quantiles", {
  toy <- toy_grid()
  bands <- default_bands(toy$coords, toy$values, n = 10)

  expect_type(bands, "double")
  expect_length(bands, 9L)
  expect_true(all(bands > 0))
  expect_equal(bands, default_bands(toy$coords, toy$values, n = 10))
})

test_that("spseg accepts smoothing as a named list", {
  data(segdata)
  geom <- sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10),
                                    c(0, 10), c(0, 0)))) |>
    sf::st_make_grid(cellsize = 1)
  grid_sf <- sf::st_sf(geometry = geom)

  result <- spseg(
    grid_sf,
    data = segdata[, 1:2],
    smoothing = list(smoothing = "equal", nrow = 10, ncol = 10),
    maxdist = 2
  )

  expect_s4_class(result, "SegSpatial")
  expect_equal(nrow(result@env), nrow(result@coords))
})

test_that("spseg accommodates deprecated top-level smoothing arguments", {
  data(segdata)
  geom <- sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10),
                                    c(0, 10), c(0, 0)))) |>
    sf::st_make_grid(cellsize = 1)
  grid_sf <- sf::st_sf(geometry = geom)

  expect_warning(
    old <- spseg(grid_sf, data = segdata[, 1:2], smoothing = "equal",
                 nrow = 10, ncol = 10, maxdist = 2),
    "Deprecated smoothing"
  )
  new <- spseg(grid_sf, data = segdata[, 1:2],
               smoothing = list(smoothing = "equal", nrow = 10, ncol = 10),
               maxdist = 2)

  expect_equal(old@env, new@env, tolerance = 1e-12)
})

test_that("SegLocal objects support spplot examples", {
  toy <- toy_grid(cols = 5:6)
  env <- localenv(toy$coords, toy$values, weighting = "biweight", power = 1)

  expect_s3_class(as(env, "sf"), "sf")
  expect_equal(nrow(as(env, "sf")), nrow(toy$coords))
  expect_s3_class(spplot(env, main = "Biweight with p = 1"), "trellis")
})

test_that("SegLocal sf coercion preserves source geometry when available", {
  geom <- sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2))),
    n = c(2, 2)
  )
  values <- matrix(c(10, 0,
                     0, 10,
                     5, 5,
                     3, 7),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("a", "b")))
  x <- sf::st_sf(a = values[, 1], b = values[, 2], geometry = geom)

  env_poly <- localenv(x, maxdist = 0)
  env_poly_sf <- as(env_poly, "sf")
  env_point <- localenv(sf::st_coordinates(sf::st_centroid(geom)), values,
                        maxdist = 0)
  env_point_sf <- as(env_point, "sf")

  expect_s3_class(env_poly_sf, "sf")
  expect_equal(as.character(unique(sf::st_geometry_type(env_poly_sf))),
               "POLYGON")
  expect_equal(nrow(env_poly_sf), nrow(x))
  expect_equal(sf::st_crs(env_poly_sf), sf::st_crs(x))
  expect_s3_class(env_point_sf, "sf")
  expect_equal(as.character(unique(sf::st_geometry_type(env_point_sf))),
               "POINT")
})

test_that("as.data.frame returns centroid coordinates and local environments", {
  toy <- toy_grid(cols = 5:6)
  env <- localenv(toy$coords, toy$values, maxdist = 0)
  env_df <- as.data.frame(env)

  expect_equal(names(env_df), c("x", "y", colnames(toy$values)))
  expect_equal(nrow(env_df), nrow(toy$coords))
  expect_equal(as.matrix(env_df[, c("x", "y")]), toy$coords,
               ignore_attr = TRUE)
  expect_equal(as.matrix(env_df[, colnames(toy$values)]), toy$values,
               ignore_attr = TRUE)
})

test_that("st_as_sf compiles selected SegResult local environments", {
  geom <- sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2))),
    n = c(2, 2)
  )
  values <- matrix(c(10, 0,
                     0, 10,
                     5, 5,
                     3, 7),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("a", "b")))
  x <- sf::st_sf(a = values[, 1], b = values[, 2], geometry = geom)
  result <- spseg(x, bands = c(0, 2), output = "localenv")

  wide <- sf::st_as_sf(result, bands = c(0, 2), columns = "a")
  long <- sf::st_as_sf(result, bands = c(0, 2), columns = "a",
                       format = "long")

  expect_s3_class(wide, "sf")
  expect_equal(names(sf::st_drop_geometry(wide)), c("bw_0_a", "bw_2_a"))
  expect_equal(nrow(wide), nrow(x))
  expect_equal(as.character(unique(sf::st_geometry_type(wide))), "POLYGON")
  expect_s3_class(long, "sf")
  expect_equal(nrow(long), nrow(x) * 2)
  expect_equal(sort(unique(long$band)), c(0, 2))
  expect_equal(names(sf::st_drop_geometry(long)), c("a", "band"))
})

test_that("as.data.frame compiles selected SegResult local environments", {
  toy <- toy_grid(cols = 5:6)
  result <- spseg(toy$coords, toy$values, bands = c(0, 2),
                  output = "localenv")

  wide <- as.data.frame(result, what = "env", bands = c(0, 2),
                        columns = "C1")
  long <- as.data.frame(result, what = "env", bands = c(0, 2),
                        columns = "C1", format = "long")

  expect_equal(names(wide), c("x", "y", "bw_0_C1", "bw_2_C1"))
  expect_equal(nrow(wide), nrow(toy$coords))
  expect_equal(names(long), c("x", "y", "band", "C1"))
  expect_equal(nrow(long), nrow(toy$coords) * 2)
  expect_equal(sort(unique(long$band)), c(0, 2))
})

test_that("test_data spatial fixtures load and feed localenv", {
  fixtures <- new.env(parent = emptyenv())
  load(testthat::test_path("fixtures", "SegAll.RData"), envir = fixtures)

  for (name in paste0("SegP", 1:5)) {
    fixture <- fixtures[[name]]
    expect_s4_class(fixture, "SpatialPolygonsDataFrame")
    expect_equal(names(fixture@data), c("grp1", "grp2"))
    expect_equal(dissim(data = fixture@data)$d, 1)

    env <- localenv(sf::st_as_sf(fixture))
    expect_s4_class(env, "SegLocal")
    expect_equal(dim(env@data), c(nrow(fixture@data), 2L))
    expect_false(anyNA(env@env))
  }
})

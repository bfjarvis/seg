test_that("localenv accepts coordinate matrices and preserves dimensions", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  env <- localenv(coords, values, useExp = FALSE, power = 0, scale = FALSE)

  expect_s4_class(env, "SegLocal")
  expect_equal(dim(env@coords), c(100L, 2L))
  expect_equal(dim(env@data), dim(values))
  expect_equal(dim(env@env), dim(values))
  expected <- matrix(colMeans(values), nrow = 100, ncol = 2, byrow = TRUE)
  dimnames(expected) <- dimnames(env@env)
  expect_equal(env@env, expected)
})

test_that("localenv default weights are invariant to coordinate units", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  env_m <- localenv(coords, values)
  env_km <- localenv(coords * 1000, values)

  expect_equal(env_km@env, env_m@env, tolerance = 1e-12)
})

test_that("localenv default bandwidth is geometric mean distance", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])
  bandwidth <- exp(mean(log(as.numeric(dist(coords)))))

  env_default <- localenv(coords, values)
  env_explicit <- localenv(coords, values, maxdist = bandwidth)

  expect_equal(env_default@env, env_explicit@env, tolerance = 1e-12)
})

test_that("automatic bandwidth is deterministic and preserves RNG state", {
  coords <- cbind(seq_len(150), sin(seq_len(150)))

  set.seed(42)
  old_seed <- .Random.seed
  bandwidth_1 <- seg:::.geometric_mean_distance(coords)
  seed_after <- .Random.seed
  bandwidth_2 <- seg:::.geometric_mean_distance(coords)

  expect_equal(seed_after, old_seed)
  expect_equal(bandwidth_2, bandwidth_1)
  expect_true(is.finite(bandwidth_1))
  expect_gt(bandwidth_1, 0)
})

test_that("localenv distance-matrix path uses normalized distances by default", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  env_coords <- localenv(coords, values)
  env_dist <- localenv(coords, values, sprel = dist(coords))

  expect_equal(env_dist@env, env_coords@env, tolerance = 1e-12)
})

test_that("localenv and spatseg vary across bandwidths", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  env_small <- localenv(coords, values, maxdist = 1)
  env_large <- localenv(coords, values, maxdist = 5)
  seg_small <- as.list(spatseg(env_small))
  seg_large <- as.list(spatseg(env_large))

  expect_false(isTRUE(all.equal(env_small@env, env_large@env)))
  expect_false(isTRUE(all.equal(seg_small$d, seg_large$d)))
})

test_that("spatseg selected methods match all-method C results", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])
  env <- localenv(coords, values, maxdist = 2)

  all_methods <- as.list(spatseg(env, useC = TRUE))
  exposure <- as.list(spatseg(env, method = "exposure", useC = TRUE))
  information <- as.list(spatseg(env, method = "information", useC = TRUE))
  diversity <- as.list(spatseg(env, method = "diversity", useC = TRUE))
  dissimilarity <- as.list(spatseg(env, method = "dissimilarity", useC = TRUE))

  expect_equal(exposure$p, all_methods$p, tolerance = 1e-12)
  expect_equal(information$h, all_methods$h, tolerance = 1e-12)
  expect_equal(diversity$r, all_methods$r, tolerance = 1e-12)
  expect_equal(dissimilarity$d, all_methods$d, tolerance = 1e-12)
})

test_that("spseg wrapper returns selected measures only", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  result <- spseg(coords, values, method = c("information", "dissimilarity"),
                  maxdist = 2)

  expect_s4_class(result, "SegSpatial")
  expect_length(result@h, 1L)
  expect_length(result@d, 1L)
  expect_length(result@r, 0L)
  expect_equal(dim(result@p), c(0L, 0L))
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

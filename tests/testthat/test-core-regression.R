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

test_that("localenv and spseg vary across bandwidths", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])

  env_small <- localenv(coords, values, maxdist = 1)
  env_large <- localenv(coords, values, maxdist = 5)
  seg_small <- as.list(spseg(env_small))
  seg_large <- as.list(spseg(env_large))

  expect_false(isTRUE(all.equal(env_small@env, env_large@env)))
  expect_false(isTRUE(all.equal(seg_small$d, seg_large$d)))
})

test_that("spseg selected methods match all-method C results", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])
  env <- localenv(coords, values, maxdist = 2)

  all_methods <- as.list(spseg(env, useC = TRUE))
  exposure <- as.list(spseg(env, method = "exposure", useC = TRUE))
  information <- as.list(spseg(env, method = "information", useC = TRUE))
  diversity <- as.list(spseg(env, method = "diversity", useC = TRUE))
  dissimilarity <- as.list(spseg(env, method = "dissimilarity", useC = TRUE))

  expect_equal(exposure$p, all_methods$p, tolerance = 1e-12)
  expect_equal(information$h, all_methods$h, tolerance = 1e-12)
  expect_equal(diversity$r, all_methods$r, tolerance = 1e-12)
  expect_equal(dissimilarity$d, all_methods$d, tolerance = 1e-12)
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

  c_result <- as.list(spseg(env, method = "information", useC = TRUE))
  r_result <- as.list(spseg(env, method = "information", useC = FALSE))

  expect_true(is.finite(c_result$h))
  expect_equal(c_result$h, r_result$h, tolerance = 1e-12)
  expect_equal(c_result$h, 1, tolerance = 1e-12)
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
    with_tol <- spseg(env, method = "information", tol = 1e-3),
    "ignored"
  )
  without_tol <- spseg(env, method = "information")

  expect_equal(as.list(with_tol)$h, as.list(without_tol)$h)
})

test_that("spatseg remains as a deprecated compatibility wrapper", {
  data(segdata)
  coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
  values <- as.matrix(segdata[, 1:2])
  env <- localenv(coords, values, maxdist = 2)

  expect_warning(
    result <- spatseg(env, method = "dissimilarity", useC = TRUE),
    "deprecated"
  )
  expect_equal(as.list(result)$d,
               as.list(spseg(env, method = "dissimilarity", useC = TRUE))$d,
               tolerance = 1e-12)
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

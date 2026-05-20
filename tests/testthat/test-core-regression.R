toy_grid <- function(cols = 1:2) {
  data(segdata)
  list(
    coords = as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5))),
    values = as.matrix(segdata[, cols])
  )
}

surface_totals_by_polygon <- function(result, source_sf) {
  zone_ids <- seg:::spseg_surface_polygon_ids(result$coords, source_sf)
  totals <- rowsum(result$data, zone_ids)
  totals[as.character(seq_len(nrow(source_sf))), , drop = FALSE]
}

test_that("spseg returns S3 index results by default", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values)

  expect_s3_class(result, "seg_result")
  expect_null(result$coords)
  expect_null(result$data)
  expect_null(result$env)
  expect_equal(result$bands, default_bands(toy$coords, toy$values))
  expect_length(result$indices$overall$d, 5L)
})

test_that("spseg can store local environments without indices", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values, bands = c(0, 2),
                  output = "localenv")

  expect_s3_class(result, "seg_result")
  expect_equal(result$bands, c(0, 2))
  expect_length(result$env, 2L)
  expect_equal(result$env[[1]], toy$values)
  expect_equal(dim(result$env[[2]]), dim(toy$values))
  expect_equal(length(result$indices), 0L)
})

test_that("spseg full output stores local environments and indices", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values, bands = c(1, 2),
                  output = "full")
  env_only <- spseg(toy$coords, toy$values, bands = 2,
                    output = "localenv")
  indices_only <- spseg(toy$coords, toy$values, bands = 2,
                        output = "indices")

  expect_s3_class(result, "seg_result")
  expect_equal(result$env[[2]], env_only$env[[1]], tolerance = 1e-12)
  expect_equal(result$indices$overall$d[["2"]],
               indices_only$indices$overall$d[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$overall$r[["2"]],
               indices_only$indices$overall$r[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$overall$h[["2"]],
               indices_only$indices$overall$h[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$overall$p[["2"]],
               indices_only$indices$overall$p[["2"]], tolerance = 1e-12)
})

test_that("indices from stored local environments match full spseg output", {
  toy <- toy_grid(cols = 1:3)
  result <- spseg(toy$coords, toy$values, bands = 2,
                  comparison = "both", output = "full")
  indices <- seg:::spseg_indices_from_env(
    data = result$data,
    env = result$env[[1]],
    measures = "all",
    comparison = "both"
  )

  expect_equal(indices$overall$d[[1]],
               result$indices$overall$d[["2"]], tolerance = 1e-12)
  expect_equal(indices$overall$r[[1]],
               result$indices$overall$r[["2"]], tolerance = 1e-12)
  expect_equal(indices$overall$h[[1]],
               result$indices$overall$h[["2"]], tolerance = 1e-12)
  expect_equal(indices$overall$p[[1]],
               result$indices$overall$p[["2"]], tolerance = 1e-12)
  expect_equal(indices$pairwise$d[[1]],
               result$indices$pairwise$d[["2"]], tolerance = 1e-12)
  expect_equal(indices$pairwise$r[[1]],
               result$indices$pairwise$r[["2"]], tolerance = 1e-12)
  expect_equal(indices$pairwise$h[[1]],
               result$indices$pairwise$h[["2"]], tolerance = 1e-12)
})

test_that("weighting schemes behave as expected for coordinates", {
  coords <- cbind(x = c(0, 1, 3), y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     10, 0),
                   ncol = 2, byrow = TRUE)
  colnames(values) <- c("a", "b")

  unweighted <- spseg(coords, values, weighting = "unweighted", bands = 1,
                      output = "localenv")
  expect_equal(unweighted$env[[1]][1, ], c(a = 5, b = 5))
  expect_equal(unweighted$env[[1]][3, ], c(a = 10, b = 0))

  inverse <- spseg(coords, values, weighting = "inverse", bands = 2,
                   power = 2, output = "localenv")
  exponential <- spseg(coords, values, weighting = "exponential", bands = 2,
                       output = "localenv")
  expect_false(isTRUE(all.equal(inverse$env, exponential$env)))
})

test_that("selected measures control stored indices", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values,
                  measures = c("information", "dissimilarity"),
                  bands = c(1, 2),
                  output = "indices")

  expect_length(result$indices$overall$h, 2L)
  expect_length(result$indices$overall$d, 2L)
  expect_length(result$indices$overall$r, 0L)
  expect_length(result$indices$overall$p, 0L)
  expect_equal(as.data.frame(result)$measure, rep(c("d", "h"), each = 2))
})

test_that("zero bandwidth produces unsmoothed indices", {
  coords <- cbind(x = 1:4, y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     5, 0,
                     0, 5),
                   ncol = 2, byrow = TRUE)
  colnames(values) <- c("a", "b")

  result <- spseg(coords, values, bands = 0, output = "full",
                  measures = "information")

  expect_equal(result$env[[1]], values)
  expect_true(is.finite(result$indices$overall$h[[1]]))
  expect_equal(result$indices$overall$h[[1]], 1, tolerance = 1e-12)
})

test_that("pairwise scalar indices agree with overall indices for two groups", {
  toy <- toy_grid(cols = 5:6)
  result <- spseg(toy$coords, toy$values, bands = 2,
                  output = "indices", comparison = "both")

  for (band in names(result$indices$overall$d)) {
    expect_equal(result$indices$pairwise$d[[band]][1, 2],
                 result$indices$overall$d[[band]], tolerance = 1e-12)
    expect_equal(result$indices$pairwise$r[[band]][1, 2],
                 result$indices$overall$r[[band]], tolerance = 1e-12)
    expect_equal(result$indices$pairwise$h[[band]][1, 2],
                 result$indices$overall$h[[band]], tolerance = 1e-12)
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

  expect_equal(result$indices$pairwise$d[["0"]]["a", "b"], 1,
               tolerance = 1e-12)
  expect_true(is.na(result$indices$pairwise$d[["0"]]["a", "empty"]))
  expect_true(is.finite(result$indices$pairwise$r[["0"]]["a", "b"]))
  expect_true(is.finite(result$indices$pairwise$h[["0"]]["a", "b"]))
  expect_equal(nrow(as.data.frame(result, what = "pairwise")), 18L)
})

test_that("sf inputs preserve geometry in local environment outputs", {
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
                       shape = "long")

  expect_s3_class(wide, "sf")
  expect_equal(names(sf::st_drop_geometry(wide)), c("bw_0_a", "bw_2_a"))
  expect_equal(as.character(unique(sf::st_geometry_type(wide))), "POLYGON")
  expect_s3_class(long, "sf")
  expect_equal(nrow(long), nrow(x) * 2)
  expect_equal(sort(unique(long$band)), c(0, 2))
})

test_that("as.data.frame compiles local environments", {
  toy <- toy_grid(cols = 5:6)
  result <- spseg(toy$coords, toy$values, bands = c(0, 2),
                  output = "localenv")

  wide <- as.data.frame(result, what = "env", bands = c(0, 2),
                        columns = "C1")
  long <- as.data.frame(result, what = "env", bands = c(0, 2),
                        columns = "C1", shape = "long")

  expect_equal(names(wide), c("x", "y", "bw_0_C1", "bw_2_C1"))
  expect_equal(nrow(wide), nrow(toy$coords))
  expect_equal(names(long), c("x", "y", "band", "C1"))
  expect_equal(nrow(long), nrow(toy$coords) * 2)
  expect_equal(sort(unique(long$band)), c(0, 2))
})

test_that("grid population surfaces preserve source polygon totals", {
  data(segdata)
  geom <- sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10),
                                    c(0, 10), c(0, 0)))) |>
    sf::st_make_grid(cellsize = 1)
  grid_sf <- sf::st_sf(geometry = geom)

  result <- spseg(
    grid_sf,
    data = segdata[, 1:2],
    surface = "grid",
    nrow = 10,
    ncol = 10,
    bands = 2,
    output = "localenv"
  )

  expect_s3_class(result, "seg_result")
  expect_s3_class(result$geometry, "sfc")
  expect_equal(as.character(unique(sf::st_geometry_type(result$geometry))),
               "POLYGON")
  surface_totals <- surface_totals_by_polygon(result, grid_sf)
  expect_equal(unname(surface_totals), unname(as.matrix(segdata[, 1:2])),
               tolerance = 1e-8)
})

test_that("pycnophylactic surfaces preserve source polygon totals", {
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

  result <- spseg(x, bands = 1, surface = "pycno", celldim = 0.5,
                  converge = 1, output = "localenv")

  expect_s3_class(result, "seg_result")
  expect_gt(nrow(result$data), nrow(x))
  surface_totals <- surface_totals_by_polygon(result, x)
  expect_equal(unname(surface_totals), unname(values), tolerance = 1e-8)
})

test_that("kNN uses count-based thresholds with fractional boundary units", {
  coords <- cbind(x = c(0, 1, 2), y = 0)
  values <- matrix(c(10, 0,
                     0, 10,
                     0, 10),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("a", "b")))

  result <- spseg(coords, values, bands = c(5, 15), neighbors = "knn",
                  weighting = "unweighted", normalize = FALSE,
                  output = "localenv")

  expect_equal(result$neighbors$type, "knn")
  expect_equal(result$neighbors$units, "population")
  expect_equal(result$env[[1]], values)
  expect_equal(result$env[[2]][1, ], c(a = 10, b = 5))
  expect_equal(result$env[[2]][2, ], c(a = 5, b = 10))
  expect_equal(result$env[[2]][3, ], c(a = 0, b = 15))
})

test_that("brute force search matches kd-tree search", {
  toy <- toy_grid(cols = 5:6)

  radius_tree <- spseg(toy$coords, toy$values, bands = c(1, 2),
                       search = "kdtree", output = "full")
  radius_brute <- spseg(toy$coords, toy$values, bands = c(1, 2),
                        search = "brute", output = "full")

  expect_equal(radius_brute$neighbors$engine, "brute")
  expect_equal(radius_brute$env, radius_tree$env, tolerance = 1e-12)
  expect_equal(radius_brute$indices, radius_tree$indices, tolerance = 1e-12)

  knn_tree <- spseg(toy$coords, toy$values, bands = c(10, 25),
                    neighbors = "knn", weighting = "unweighted",
                    normalize = FALSE, search = "kdtree", output = "full")
  knn_brute <- spseg(toy$coords, toy$values, bands = c(10, 25),
                     neighbors = "knn", weighting = "unweighted",
                     normalize = FALSE, search = "brute", output = "full")

  expect_equal(knn_brute$neighbors$engine, "brute")
  expect_equal(knn_brute$env, knn_tree$env, tolerance = 1e-12)
  expect_equal(knn_brute$indices, knn_tree$indices, tolerance = 1e-12)
})

toy_grid <- function(cols = 1:2) {
  data(segdata)
  list(
    coords = as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5))),
    values = as.matrix(segdata[, cols])
  )
}

test_that("spseg returns S3 index results by default", {
  toy <- toy_grid()
  result <- spseg(toy$coords, toy$values)

  expect_s3_class(result, "seg_result")
  expect_null(result$coords)
  expect_null(result$data)
  expect_null(result$env)
  expect_equal(result$bands, default_bands(toy$coords, toy$values))
  expect_length(result$indices$multigroup$d, length(result$bands))
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
  expect_equal(result$indices$multigroup$d[["2"]],
               indices_only$indices$multigroup$d[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$multigroup$r[["2"]],
               indices_only$indices$multigroup$r[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$multigroup$h[["2"]],
               indices_only$indices$multigroup$h[["2"]], tolerance = 1e-12)
  expect_equal(result$indices$multigroup$p[["2"]],
               indices_only$indices$multigroup$p[["2"]], tolerance = 1e-12)
})

test_that("indices from stored local environments match full spseg output", {
  toy <- toy_grid(cols = 1:3)
  result <- spseg(toy$coords, toy$values, bands = 2,
                  scope = "both", output = "full")
  indices <- seg:::spseg_indices_from_env(
    data = result$data,
    env = result$env[[1]],
    measures = "all",
    scope = "both"
  )

  expect_equal(indices$multigroup$d[[1]],
               result$indices$multigroup$d[["2"]], tolerance = 1e-12)
  expect_equal(indices$multigroup$r[[1]],
               result$indices$multigroup$r[["2"]], tolerance = 1e-12)
  expect_equal(indices$multigroup$h[[1]],
               result$indices$multigroup$h[["2"]], tolerance = 1e-12)
  expect_equal(indices$multigroup$p[[1]],
               result$indices$multigroup$p[["2"]], tolerance = 1e-12)
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

  expect_length(result$indices$multigroup$h, 2L)
  expect_length(result$indices$multigroup$d, 2L)
  expect_length(result$indices$multigroup$r, 0L)
  expect_length(result$indices$multigroup$p, 0L)
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
  expect_true(is.finite(result$indices$multigroup$h[[1]]))
  expect_equal(result$indices$multigroup$h[[1]], 1, tolerance = 1e-12)
})

test_that("pairwise scalar indices agree with multigroup indices for two groups", {
  toy <- toy_grid(cols = 5:6)
  result <- spseg(toy$coords, toy$values, bands = 2,
                  output = "indices", scope = "both")

  for (band in names(result$indices$multigroup$d)) {
    expect_equal(result$indices$pairwise$d[[band]][1, 2],
                 result$indices$multigroup$d[[band]], tolerance = 1e-12)
    expect_equal(result$indices$pairwise$r[[band]][1, 2],
                 result$indices$multigroup$r[[band]], tolerance = 1e-12)
    expect_equal(result$indices$pairwise$h[[band]][1, 2],
                 result$indices$multigroup$h[[band]], tolerance = 1e-12)
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
                  scope = "pairwise",
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

  checked <- suppressMessages(seg:::chksegdata(grid_sf, segdata[, 1:2]))
  surface <- seg:::spseg_surface(
    x = grid_sf,
    coords = checked$coords,
    data = checked$data,
    surface = "grid",
    args = list(nrow = 10, ncol = 10),
    verbose = FALSE
  )

  result <- spseg(grid_sf, data = segdata[, 1:2], surface = "grid",
                  nrow = 10, ncol = 10, bands = 2, output = "localenv")
  expect_s3_class(result, "seg_result")
  expect_null(result$geometry)
  surface_totals <- rowsum(surface$data, surface$id)
  expect_equal(unname(surface_totals), unname(as.matrix(segdata[, 1:2])),
               tolerance = 1e-8)

  result_with_points <- spseg(
    grid_sf,
    data = segdata[, 1:2],
    surface = "grid",
    surface_geometry = "points",
    nrow = 10,
    ncol = 10,
    bands = 2,
    output = "localenv"
  )
  expect_s3_class(result_with_points$geometry, "sfc")
  expect_equal(
    as.character(unique(sf::st_geometry_type(result_with_points$geometry))),
    "POINT"
  )

  result_with_geometry <- spseg(
    grid_sf,
    data = segdata[, 1:2],
    surface = "grid",
    surface_geometry = TRUE,
    nrow = 10,
    ncol = 10,
    bands = 2,
    output = "localenv"
  )
  expect_s3_class(result_with_geometry$geometry, "sfc")
  expect_equal(
    as.character(unique(sf::st_geometry_type(result_with_geometry$geometry))),
    "POLYGON"
  )
})

test_that("raw, grid, and pycno surfaces preserve overall group totals", {
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
  checked <- suppressMessages(seg:::chksegdata(x[, c("a", "b")]))

  raw <- seg:::spseg_surface(
    x = x[, c("a", "b")],
    coords = checked$coords,
    data = checked$data,
    surface = "raw",
    args = list(),
    verbose = FALSE
  )
  grid <- seg:::spseg_surface(
    x = x[, c("a", "b")],
    coords = checked$coords,
    data = checked$data,
    surface = "grid",
    args = list(celldim = 0.5),
    verbose = FALSE
  )
  pycno <- seg:::spseg_surface(
    x = x[, c("a", "b")],
    coords = checked$coords,
    data = checked$data,
    surface = "pycno",
    args = list(celldim = 0.5, converge = 1),
    verbose = FALSE
  )

  expected <- colSums(values)
  expect_equal(unname(colSums(raw$data)), unname(expected), tolerance = 1e-8)
  expect_equal(unname(colSums(grid$data)), unname(expected), tolerance = 1e-8)
  expect_equal(unname(colSums(pycno$data)), unname(expected), tolerance = 1e-8)
})

test_that("surface audit retains grid and pycno intermediates", {
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

  minimal <- spseg(
    x,
    bands = 1,
    surface = "grid",
    nrow = 4,
    ncol = 4,
    output = "indices"
  )
  audited <- spseg(
    x,
    bands = 1,
    surface = "pycno",
    surface_audit = TRUE,
    nrow = 4,
    ncol = 4,
    converge = 1,
    output = "indices"
  )

  expect_null(minimal$surface_info$audit)
  expect_equal(minimal$surface_info$fallback$count, 0L)
  expect_null(audited$data)
  expect_equal(audited$surface_info$type, "pycno")
  expect_equal(audited$surface_info$grid$retained_cells, 16L)
  expect_equal(
    unname(audited$surface_info$audit$source$data),
    unname(values)
  )
  expect_s3_class(audited$surface_info$audit$source$geometry, "sfc")
  expect_equal(nrow(audited$surface_info$audit$cells), 16L)
  expect_equal(dim(audited$surface_info$audit$initial_counts), c(16L, 2L))
  expect_equal(dim(audited$surface_info$audit$smoothed_counts), c(16L, 2L))
  expect_equal(
    colSums(audited$surface_info$audit$initial_counts),
    colSums(audited$surface_info$audit$smoothed_counts),
    tolerance = 1e-8
  )
})

test_that("grid surfaces assign missing non-empty polygons to fallback cells", {
  tiny <- sf::st_polygon(list(rbind(c(0, 0), c(0.1, 0), c(0.1, 0.1),
                                    c(0, 0.1), c(0, 0))))
  large <- sf::st_polygon(list(rbind(c(4, 4), c(10, 4), c(10, 10),
                                     c(4, 10), c(4, 4))))
  x <- sf::st_sf(
    a = c(10, 20),
    b = c(0, 5),
    geometry = sf::st_sfc(tiny, large)
  )
  checked <- suppressMessages(seg:::chksegdata(x[, c("a", "b")]))

  expect_warning(
    result <- seg:::spseg_surface(
      x = x[, c("a", "b")],
      coords = checked$coords,
      data = checked$data,
      surface = "grid",
      args = list(nrow = 2, ncol = 2),
      audit = TRUE,
      verbose = FALSE
    ),
    "1 populated source polygon"
  )

  expect_true(1L %in% result$fallback$ids)
  expect_equal(result$surface_info$fallback$count, 1L)
  expect_equal(result$surface_info$audit$fallback$source_id, 1L)
  expect_equal(
    as.numeric(result$surface_info$audit$fallback[1, c("a", "b")]),
    c(10, 0)
  )
  expect_equal(
    colSums(result$data),
    colSums(as.matrix(sf::st_drop_geometry(x[, c("a", "b")]))),
    tolerance = 1e-8
  )
  assigned <- !is.na(result$id)
  surface_totals <- rowsum(result$data[assigned, , drop = FALSE],
                           result$id[assigned])
  expect_equal(
    unname(surface_totals["2", , drop = FALSE]),
    unname(as.matrix(sf::st_drop_geometry(x[2, c("a", "b")]))),
    tolerance = 1e-8
  )

  expect_warning(
    pycno <- seg:::spseg_surface(
      x = x[, c("a", "b")],
      coords = checked$coords,
      data = checked$data,
      surface = "pycno",
      args = list(nrow = 2, ncol = 2, converge = 1),
      verbose = FALSE
    ),
    "not preserved as separate zonal constraints"
  )
  expect_equal(
    unname(colSums(pycno$data)),
    unname(colSums(as.matrix(sf::st_drop_geometry(x[, c("a", "b")])))),
    tolerance = 1e-8
  )
})

test_that("grid fallback preserves the locations of missing polygon counts", {
  tiny_1 <- sf::st_polygon(list(rbind(c(0, 0), c(0.2, 0), c(0.2, 1),
                                      c(0, 1), c(0, 0))))
  tiny_2 <- sf::st_polygon(list(rbind(c(0, 1), c(0.2, 1), c(0.2, 2),
                                      c(0, 2), c(0, 1))))
  large_1 <- sf::st_polygon(list(rbind(c(0.2, 0), c(2, 0), c(2, 2),
                                       c(0.2, 2), c(0.2, 0))))
  large_2 <- sf::st_polygon(list(rbind(c(2, 0), c(4, 0), c(4, 2),
                                       c(2, 2), c(2, 0))))
  x <- sf::st_sf(
    a = c(10, 0, 20, 0),
    b = c(0, 6, 4, 30),
    geometry = sf::st_sfc(tiny_1, tiny_2, large_1, large_2)
  )
  checked <- suppressMessages(seg:::chksegdata(x[, c("a", "b")]))

  expect_warning(
    result <- seg:::spseg_surface(
      x = x[, c("a", "b")],
      coords = checked$coords,
      data = checked$data,
      surface = "grid",
      args = list(nrow = 2, ncol = 4),
      verbose = FALSE
    ),
    "2 populated source polygon"
  )

  expect_true(all(c(1L, 2L) %in% result$fallback$ids))
  fallback_cells <- result$fallback$cells[
    match(c(1L, 2L), result$fallback$ids)
  ]
  expect_length(unique(fallback_cells), 2L)
  expect_equal(
    unname(colSums(result$data)),
    unname(colSums(as.matrix(sf::st_drop_geometry(x[, c("a", "b")])))),
    tolerance = 1e-8
  )

  base_large_1 <- c(a = 5, b = 1)
  expect_equal(
    unname(result$data[fallback_cells[1], ]),
    unname(base_large_1 + c(a = 10, b = 0)),
    tolerance = 1e-8
  )
  expect_equal(
    unname(result$data[fallback_cells[2], ]),
    unname(base_large_1 + c(a = 0, b = 6)),
    tolerance = 1e-8
  )
  other_large_1_cells <- setdiff(which(result$id == 3L), fallback_cells)
  expect_equal(
    unname(result$data[other_large_1_cells, , drop = FALSE]),
    unname(matrix(base_large_1, nrow = length(other_large_1_cells),
                  ncol = 2, byrow = TRUE)),
    tolerance = 1e-8
  )

  expect_warning(
    pycno <- seg:::spseg_surface(
      x = x[, c("a", "b")],
      coords = checked$coords,
      data = checked$data,
      surface = "pycno",
      args = list(nrow = 2, ncol = 4, converge = 1),
      verbose = FALSE
    ),
    "not preserved as separate zonal constraints"
  )
  expect_equal(
    unname(rowsum(pycno$data, pycno$id)),
    unname(rowsum(result$data, result$id)),
    tolerance = 1e-8
  )
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

  checked <- suppressMessages(seg:::chksegdata(x[, c("a", "b")]))
  surface <- seg:::spseg_surface(
    x = x[, c("a", "b")],
    coords = checked$coords,
    data = checked$data,
    surface = "pycno",
    args = list(celldim = 0.5, converge = 1),
    verbose = FALSE
  )
  result <- spseg(x, bands = 1, surface = "pycno", celldim = 0.5,
                  converge = 1, output = "localenv")

  expect_s3_class(result, "seg_result")
  expect_gt(nrow(result$data), nrow(x))
  surface_totals <- rowsum(surface$data, surface$id)
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

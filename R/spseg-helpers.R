# ------------------------------------------------------------------------------
# Internal helpers for spseg()
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg_measures <- function(measures) {
  measures <- match.arg(measures, c("exposure", "information", "diversity",
                                    "dissimilarity", "all"), several.ok = TRUE)
  if ("all" %in% measures)
    measures <- c("exposure", "information", "diversity", "dissimilarity")

  measures
}

spseg_empty_results <- function() {
  list(p = matrix(0, nrow = 0, ncol = 0), h = numeric(), r = numeric(),
       d = numeric())
}

spseg_measure_flags <- function(measures) {
  c("exposure" %in% measures, "information" %in% measures,
    "diversity" %in% measures, "dissimilarity" %in% measures)
}

.geometric_mean_distance <- function(x, max_points = 500) {
  n <- nrow(x)
  if (n > max_points) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    set.seed(1)
    on.exit({
      if (has_seed)
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    x <- x[sample.int(n, max_points), , drop = FALSE]
  }

  d <- as.numeric(dist(x))
  d <- d[d > 0]
  if (length(d) == 0)
    return(0)

  exp(mean(log(d)))
}

default_bands <- function(x, data = NULL, n = 500,
                          probs = seq(0.1, 0.5, 0.1),
                          weighted = TRUE) {
  tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords
  data <- tmp$data
  n_points <- nrow(coords)
  sample_n <- min(n_points, n)

  weights <- NULL
  if (weighted) {
    weights <- rowSums(data)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0)
      weights <- NULL
  }

  if (n_points > sample_n) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    set.seed(1)
    on.exit({
      if (has_seed)
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    coords <- coords[sample.int(n_points, sample_n, prob = weights), ,
                     drop = FALSE]
  }

  d <- as.numeric(dist(coords))
  d <- d[d > 0]
  if (length(d) == 0)
    stop("failed to calculate default bands from the input coordinates",
         call. = FALSE)

  as.numeric(quantile(d, probs = probs, names = FALSE, type = 8))
}

.restore_env_dimnames <- function(env, data) {
  colnames(env) <- colnames(data)
  rownames(env) <- rownames(data)
  env
}

spseg_weighting_id <- function(weighting) {
  match(weighting, c("unweighted", "biweight", "inverse", "exponential")) - 1L
}

spseg_dots <- function(...) {
  dots <- list(...)
  if ("tol" %in% names(dots)) {
    message("'tol' is ignored by spseg(); zero entropy terms are handled directly")
    dots$tol <- NULL
  }

  dots
}

spseg_output <- function(output) {
  match.arg(output[1], c("indices", "full", "localenv", "legacy"))
}

spseg_comparison <- function(comparison) {
  match.arg(comparison, c("overall", "pairwise", "both"))
}

spseg_comparison_flags <- function(comparison) {
  comparison <- spseg_comparison(comparison)
  c(overall = comparison %in% c("overall", "both"),
    pairwise = comparison %in% c("pairwise", "both"))
}

spseg_neighbors <- function(neighbors) {
  match.arg(neighbors, c("radius", "knn"))
}

spseg_neighbor_id <- function(neighbors) {
  match(spseg_neighbors(neighbors), c("radius", "knn")) - 1L
}

spseg_search <- function(search) {
  match.arg(search, c("kdtree", "brute"))
}

spseg_search_id <- function(search) {
  match(spseg_search(search), c("kdtree", "brute")) - 1L
}

spseg_bands <- function(maxdist, bands, output, x, data) {
  if (!is.null(bands) && !is.null(maxdist))
    stop("supply only one of 'bands' or 'maxdist'", call. = FALSE)
  if (!is.null(bands))
    return(as.numeric(bands))
  if (!is.null(maxdist))
    return(as.numeric(maxdist))
  if (identical(output, "legacy"))
    return(NULL)

  default_bands(x, data = data)
}

seg_engine_coords <- function(coords, data, bands, power, weighting, normalize,
                              measures, comparison = "overall",
                              neighbors = "radius", search = "kdtree",
                              keep_env,
                              keep_indices) {
  if (nrow(coords) != nrow(data))
    stop("'data' must have the same number of rows as 'coords'", call. = FALSE)

  xval <- coords[, 1]
  yval <- coords[, 2]
  measures <- if (keep_indices) spseg_measures(measures) else character()
  measure_flags <- spseg_measure_flags(measures)
  comparison_flags <- spseg_comparison_flags(comparison)
  neighbors_id <- spseg_neighbor_id(neighbors)
  search_id <- spseg_search_id(search)
  out <- seg_engine_cpp(
    x = xval,
    y = yval,
    data = data,
    bands = bands,
    power = power,
    weighting = as.integer(spseg_weighting_id(weighting)),
    normalize = as.integer(normalize),
    measures = as.integer(measure_flags),
    comparison = as.integer(comparison_flags),
    keep_env = as.integer(keep_env),
    keep_indices = as.integer(keep_indices),
    neighbors = as.integer(neighbors_id),
    search = as.integer(search_id)
  )

  if (!is.null(out$env)) {
    out$env <- lapply(out$env, .restore_env_dimnames, data = data)
  }
  if (!is.null(out$indices$overall$p)) {
    out$indices$overall$p <- lapply(out$indices$overall$p, function(p) {
      rownames(p) <- colnames(p) <- colnames(data)
      p
    })
  }
  for (nm in intersect(c("d", "r", "h"), names(out$indices$pairwise))) {
    out$indices$pairwise[[nm]] <- lapply(out$indices$pairwise[[nm]], function(x) {
      rownames(x) <- colnames(x) <- colnames(data)
      x
    })
  }

  out
}

spseg_indices_from_engine <- function(indices, bands) {
  if (length(indices) == 0)
    return(list())
  for (group in intersect(c("overall", "pairwise"), names(indices))) {
    if (length(indices[[group]]$d) > 0)
      indices[[group]]$d <- setNames(indices[[group]]$d, bands)
    if (length(indices[[group]]$r) > 0)
      indices[[group]]$r <- setNames(indices[[group]]$r, bands)
    if (length(indices[[group]]$h) > 0)
      indices[[group]]$h <- setNames(indices[[group]]$h, bands)
    if (length(indices[[group]]$p) > 0)
      indices[[group]]$p <- setNames(indices[[group]]$p, bands)
  }
  indices
}

spseg_legacy_from_result <- function(result) {
  indices <- result@indices$overall
  p <- if (length(indices$p) > 0) indices$p[[1]] else
    matrix(0, nrow = 0, ncol = 0)
  env <- if (!is.null(result@env)) result@env[[1]] else result@data

  SegSpatial(indices$d, indices$r, indices$h, p,
             result@coords, result@data, env, result@proj4string)
}

spseg_result <- function(coords, data, env, bands, indices, measures,
                         comparison, weighting, power, normalize, proj4string,
                         output, call, geometry = NULL,
                         neighbors = "radius", search = "kdtree") {
  keep_inputs <- !identical(output, "indices")
  neighbors <- spseg_neighbors(neighbors)
  search <- spseg_search(search)
  SegResult(
    coords = if (keep_inputs) coords else NULL,
    data = if (keep_inputs) data else NULL,
    env = env,
    geometry = if (keep_inputs) geometry else NULL,
    bands = bands,
    indices = indices,
    measures = measures,
    weighting = weighting,
    power = power,
    normalize = normalize,
    neighbors = list(type = neighbors, values = bands,
                     units = if (neighbors == "knn") "population" else
                       "distance",
                     engine = search, comparison = comparison),
    proj4string = proj4string,
    call = call
  )
}

spseg_prepare_localenv <- function(env, negative.rm) {
  dd <- env@data
  ee <- env@env

  negIDs <- rowSums(dd < 0) > 0
  if (sum(negIDs) > 0) {
    if (negative.rm) {
      warning("Rows with negative values have been removed", call. = FALSE)
      dd <- dd[!negIDs, ]
      ee <- ee[!negIDs, ]
    } else {
      warning("Negative values replaced with zero", call. = FALSE)
      dd[dd < 0] <- 0
      ee[ee < 0] <- 0
    }
  }

  negIDs <- rowSums(ee < 0) > 0
  if (sum(negIDs) > 0) {
    if (negative.rm) {
      warning("Rows with negative values removed", call. = FALSE)
      dd <- dd[!negIDs, ]
      ee <- ee[!negIDs, ]
    } else {
      warning("Negative values replaced with zero", call. = FALSE)
      ee[ee < 0] <- 0
      dd[dd < 0] <- 0
    }
  }

  list(data = dd, env = ee)
}

spseg_restore_index_dimnames <- function(indices, data) {
  if (length(indices$overall$p) > 0) {
    indices$overall$p <- lapply(indices$overall$p, function(p) {
      rownames(p) <- colnames(p) <- colnames(data)
      p
    })
  }
  for (nm in intersect(c("d", "r", "h"), names(indices$pairwise))) {
    indices$pairwise[[nm]] <- lapply(indices$pairwise[[nm]], function(x) {
      rownames(x) <- colnames(x) <- colnames(data)
      x
    })
  }
  indices
}

spseg_indices_from_env <- function(data, env, measures, comparison) {
  measures <- spseg_measures(measures)
  indices <- seg_indices_env_cpp(
    data = data,
    env = env,
    measures = as.integer(spseg_measure_flags(measures)),
    comparison = as.integer(spseg_comparison_flags(comparison))
  )
  spseg_restore_index_dimnames(indices, data)
}

spseg_indices_from_localenv <- function(env, measures, useC, negative.rm,
                                        comparison = "overall") {
  measures <- spseg_measures(measures)
  prepared <- spseg_prepare_localenv(env, negative.rm)
  spseg_indices_from_env(prepared$data, prepared$env, measures, comparison)
}

spseg_from_localenv <- function(env, measures, useC, negative.rm) {
  indices <- spseg_indices_from_localenv(env, measures, useC, negative.rm,
                                         comparison = "overall")
  results <- indices$overall
  p <- if (length(results$p) > 0) results$p[[1]] else
    matrix(0, nrow = 0, ncol = 0)

  SegSpatial(results$d, results$r, results$h, p,
             env@coords, env@data, env@env, env@proj4string)
}

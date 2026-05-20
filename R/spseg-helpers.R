# ------------------------------------------------------------------------------
# Internal helpers for spseg()
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg_measures <- function(measures) {
  measures <- match.arg(
    measures,
    c("exposure", "information", "diversity", "dissimilarity", "all"),
    several.ok = TRUE
  )
  if ("all" %in% measures) {
    measures <- c("exposure", "information", "diversity", "dissimilarity")
  }

  measures
}

spseg_measure_flags <- function(measures) {
  c(
    "exposure" %in% measures,
    "information" %in% measures,
    "diversity" %in% measures,
    "dissimilarity" %in% measures
  )
}

#' Default Distance Bands
#'
#' Calculates default distance bands from sampled pairwise distances.
#'
#' This helper samples up to `n` points, calculates pairwise distances among
#' the sampled points, and returns the requested quantiles. By default,
#' sampled points are selected with probability proportional to their total
#' population, and only the lower half of the distance distribution is used.
#' It is used by [spseg()] when no bandwidths are supplied.
#'
#' @param x A numeric coordinate matrix, data frame of coordinates, or
#'   [sf::sf] object accepted by [spseg()].
#' @param data Optional population data. If omitted and `x` contains attached
#'   data, those data are used.
#' @param n Maximum number of points to sample.
#' @param probs Probabilities passed to [stats::quantile()]. The default
#'   returns the first through fifth deciles.
#' @param weighted Logical. If TRUE, sample points using row sums of `data`
#'   as population weights.
#'
#' @return A numeric vector of distance bands.
#'
#' @seealso [spseg()]
#'
#' @export
default_bands <- function(
  x,
  data = NULL,
  n = 500,
  probs = c(0.01,0.05,seq(0.1,0.7,0.1)),
  weighted = TRUE
) {
  tmp <- if (is.null(data)) {
    suppressMessages(chksegdata(x))
  } else {
    suppressMessages(chksegdata(x, data))
  }
  coords <- tmp$coords
  data <- tmp$data
  n_points <- nrow(coords)
  sample_n <- min(n_points, n)

  weights <- NULL
  if (weighted) {
    weights <- rowSums(data)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) {
      weights <- NULL
    }
  }

  if (n_points > sample_n) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    set.seed(1)
    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    coords <- coords[
      sample.int(n_points, sample_n, prob = weights),
      ,
      drop = FALSE
    ]
  }

  d <- as.numeric(dist(coords))
  d <- d[d > 0]
  if (length(d) == 0) {
    stop(
      "failed to calculate default bands from the input coordinates",
      call. = FALSE
    )
  }

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

spseg_scope_flags <- function(scope) {
  c(
    multigroup = scope %in% c("multigroup", "both"),
    pairwise = scope %in% c("pairwise", "both")
  )
}

seg_engine_coords <- function(
  coords,
  data,
  bands,
  power,
  weighting,
  normalize,
  measures,
  scope = "both",
  neighbors = "radius",
  search = "kdtree",
  keep_env,
  keep_indices
) {
  if (nrow(coords) != nrow(data)) {
    stop("'data' must have the same number of rows as 'coords'", call. = FALSE)
  }

  xval <- coords[, 1]
  yval <- coords[, 2]
  measures <- if (keep_indices) spseg_measures(measures) else character()
  measure_flags <- spseg_measure_flags(measures)
  scope_flags <- spseg_scope_flags(scope)
  neighbors_id <- match(neighbors, c("radius", "knn")) - 1L
  search_id <- match(search, c("kdtree", "brute")) - 1L
  out <- seg_engine_cpp(
    x = xval,
    y = yval,
    data = data,
    bands = bands,
    power = power,
    weighting = as.integer(spseg_weighting_id(weighting)),
    normalize = as.integer(normalize),
    measures = as.integer(measure_flags),
    scope = as.integer(scope_flags),
    keep_env = as.integer(keep_env),
    keep_indices = as.integer(keep_indices),
    neighbors = as.integer(neighbors_id),
    search = as.integer(search_id)
  )

  if (!is.null(out$env)) {
    out$env <- lapply(out$env, .restore_env_dimnames, data = data)
  }
  if (!is.null(out$indices$multigroup$p)) {
    out$indices$multigroup$p <- lapply(out$indices$multigroup$p, function(p) {
      rownames(p) <- colnames(p) <- colnames(data)
      p
    })
  }
  for (nm in intersect(c("d", "r", "h"), names(out$indices$pairwise))) {
    out$indices$pairwise[[nm]] <- lapply(
      out$indices$pairwise[[nm]],
      function(x) {
        rownames(x) <- colnames(x) <- colnames(data)
        x
      }
    )
  }

  out
}

spseg_indices_from_engine <- function(indices, bands) {
  if (length(indices) == 0) {
    return(list())
  }
  for (group in intersect(c("multigroup", "pairwise"), names(indices))) {
    if (length(indices[[group]]$d) > 0) {
      indices[[group]]$d <- setNames(indices[[group]]$d, bands)
    }
    if (length(indices[[group]]$r) > 0) {
      indices[[group]]$r <- setNames(indices[[group]]$r, bands)
    }
    if (length(indices[[group]]$h) > 0) {
      indices[[group]]$h <- setNames(indices[[group]]$h, bands)
    }
    if (length(indices[[group]]$p) > 0) {
      indices[[group]]$p <- setNames(indices[[group]]$p, bands)
    }
  }
  indices
}

spseg_indices_from_env <- function(data, env, measures, scope = "both") {
  measures <- spseg_measures(measures)
  indices <- seg_indices_env_cpp(
    data = data,
    env = env,
    measures = as.integer(spseg_measure_flags(measures)),
    scope = as.integer(spseg_scope_flags(scope))
  )
  if (length(indices$multigroup$p) > 0) {
    indices$multigroup$p <- lapply(indices$multigroup$p, function(p) {
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

# ------------------------------------------------------------------------------
# S3 result object and methods for spseg()
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

new_seg_result <- function(coords = NULL, data = NULL, env = NULL, bands,
                           indices = list(), measures = character(),
                           weighting = character(), power = numeric(),
                           normalize = logical(),
                           neighbors = list(type = "radius"),
                           geometry = NULL,
                           crs = st_crs(as.character(NA)),
                           surface = "raw",
                           output = "indices",
                           call = NULL) {
  x <- list(
    coords = coords,
    data = data,
    env = env,
    geometry = geometry,
    bands = bands,
    indices = indices,
    measures = measures,
    weighting = weighting,
    power = power,
    normalize = normalize,
    neighbors = neighbors,
    crs = st_crs(crs),
    surface = surface,
    output = output,
    call = call
  )
  validate_seg_result(x)
  class(x) <- "seg_result"
  x
}

validate_seg_result <- function(x) {
  if (!is.null(x$coords) &&
      (!is.matrix(x$coords) || ncol(x$coords) != 2 || !is.numeric(x$coords)))
    stop("'coords' must be NULL or a numeric matrix with two columns",
         call. = FALSE)
  if (!is.null(x$data) &&
      (!is.matrix(x$data) || ncol(x$data) < 2 || !is.numeric(x$data)))
    stop("'data' must be NULL or a numeric matrix with at least two columns",
         call. = FALSE)
  if (!is.null(x$coords) && !is.null(x$data) &&
      nrow(x$coords) != nrow(x$data))
    stop("'data' must have the same number of rows as 'coords'",
         call. = FALSE)
  if (!is.null(x$env)) {
    if (!is.list(x$env))
      stop("'env' must be NULL or a list of matrices", call. = FALSE)
    if (length(x$env) != length(x$bands))
      stop("'env' must have one element per band", call. = FALSE)
  }
  if (!is.null(x$geometry) &&
      (!inherits(x$geometry, "sfc") ||
       (!is.null(x$coords) && length(x$geometry) != nrow(x$coords))))
    stop("'geometry' must be NULL or an sfc object with one feature per row",
         call. = FALSE)
  if (length(x$bands) == 0 || any(!is.finite(x$bands)))
    stop("'bands' must contain finite numeric values", call. = FALSE)
  invisible(x)
}

print.seg_result <- function(x, ...) {
  cat("\n\tSpatial segregation result\n\n")
  cat("Bands                 :", length(x$bands), "\n")
  cat("Neighborhood type     :", x$neighbors$type %||% "radius", "\n")
  cat("Neighborhood units    :", x$neighbors$units %||% "distance", "\n")
  cat("Search engine         :", x$neighbors$engine %||% "kdtree", "\n")
  cat("Surface               :", x$surface %||% "raw", "\n")
  cat("Local environments    :", if (is.null(x$env)) "not stored" else "stored", "\n")
  cat("Segregation indices   :", if (length(x$indices) == 0) "not stored" else "stored", "\n")
  if (length(x$weighting) > 0) {
    cat("Weighting             :", x$weighting, "\n")
    cat("Power                 :", x$power, "\n")
    cat("Normalize distances   :", x$normalize, "\n")
  }
  invisible(x)
}

as.list.seg_result <- function(x, ...) {
  unclass(x)
}

.seg_env_subset <- function(env, columns = NULL) {
  if (is.null(columns))
    return(env)
  if (is.numeric(columns))
    return(env[, columns, drop = FALSE])
  if (!all(columns %in% colnames(env)))
    stop("'columns' must name columns in the local environment matrix",
         call. = FALSE)
  env[, columns, drop = FALSE]
}

.seg_env_sf <- function(env, coords, geometry, crs, columns = NULL) {
  out <- data.frame(.seg_env_subset(env, columns))
  if (!is.null(geometry)) {
    st_sf(out, geometry = geometry)
  } else {
    out$x <- coords[, 1]
    out$y <- coords[, 2]
    st_as_sf(out, coords = c("x", "y"), crs = st_crs(crs))
  }
}

st_as_sf.seg_result <- function(x, ..., bands = NULL, columns = NULL,
                                shape = c("wide", "long")) {
  shape <- match.arg(shape)
  if (is.null(x$env))
    stop("'x' does not contain stored local environments", call. = FALSE)
  if (is.null(x$coords))
    stop("'x' does not contain coordinates; use output = 'full' or 'localenv'",
         call. = FALSE)

  band_index <- if (is.null(bands)) {
    seq_along(x$bands)
  } else {
    match(bands, x$bands)
  }
  if (anyNA(band_index))
    stop("'bands' must match values in x$bands", call. = FALSE)

  if (shape == "long") {
    out <- do.call(rbind, lapply(band_index, function(i) {
      y <- .seg_env_sf(x$env[[i]], x$coords, x$geometry, x$crs,
                       columns = columns)
      y$band <- x$bands[i]
      y
    }))
    rownames(out) <- NULL
    return(out)
  }

  env <- do.call(cbind, lapply(band_index, function(i) {
    y <- .seg_env_subset(x$env[[i]], columns)
    colnames(y) <- paste0(make.names(paste0("bw_", x$bands[i])), "_",
                          colnames(y))
    y
  }))
  .seg_env_sf(env, x$coords, x$geometry, x$crs)
}

as.data.frame.seg_result <- function(x, row.names = NULL, optional = FALSE,
                                     ..., what = c("indices", "pairwise",
                                                   "exposure", "env"),
                                     bands = NULL, columns = NULL,
                                     shape = c("wide", "long")) {
  what <- match.arg(what)
  shape <- match.arg(shape)

  if (what == "env") {
    if (is.null(x$env))
      stop("'x' does not contain stored local environments", call. = FALSE)
    if (is.null(x$coords))
      stop("'x' does not contain coordinates; use output = 'full' or 'localenv'",
           call. = FALSE)

    band_index <- if (is.null(bands)) seq_along(x$bands) else match(bands, x$bands)
    if (anyNA(band_index))
      stop("'bands' must match values in x$bands", call. = FALSE)

    if (shape == "long") {
      out <- do.call(rbind, lapply(band_index, function(i) {
        data.frame(x = x$coords[, 1], y = x$coords[, 2], band = x$bands[i],
                   .seg_env_subset(x$env[[i]], columns),
                   check.names = !optional)
      }))
      rownames(out) <- NULL
      if (!is.null(x$geometry))
        out$geometry <- rep(x$geometry, times = length(band_index))
      if (!is.null(row.names))
        rownames(out) <- row.names
      return(out)
    }

    env <- do.call(cbind, lapply(band_index, function(i) {
      y <- .seg_env_subset(x$env[[i]], columns)
      colnames(y) <- paste0(make.names(paste0("bw_", x$bands[i])), "_",
                            colnames(y))
      y
    }))
    out <- data.frame(x = x$coords[, 1], y = x$coords[, 2], env,
                      check.names = !optional)
    if (!is.null(x$geometry))
      out$geometry <- x$geometry
    if (!is.null(row.names))
      rownames(out) <- row.names
    return(out)
  }

  if (length(x$indices) == 0)
    stop("'x' does not contain segregation indices", call. = FALSE)

  if (what == "indices") {
    overall <- x$indices$overall
    out <- data.frame()
    scalar_names <- intersect(c("d", "r", "h"), names(overall))
    for (nm in scalar_names) {
      values <- overall[[nm]]
      if (length(values) > 0)
        out <- rbind(out, data.frame(comparison = "overall",
                                     band = x$bands, measure = nm,
                                     value = values))
    }
    rownames(out) <- NULL
    return(out)
  }

  if (what == "pairwise") {
    pairwise <- x$indices$pairwise
    out <- data.frame()
    scalar_names <- intersect(c("d", "r", "h"), names(pairwise))
    for (nm in scalar_names) {
      mats <- pairwise[[nm]]
      if (length(mats) == 0)
        next
      out <- rbind(out, do.call(rbind, lapply(seq_along(mats), function(i) {
        mat <- mats[[i]]
        idx <- which(upper.tri(mat), arr.ind = TRUE)
        data.frame(comparison = "pairwise", band = x$bands[i], measure = nm,
                   group_a = rownames(mat)[idx[, 1]],
                   group_b = colnames(mat)[idx[, 2]],
                   value = mat[idx], row.names = NULL)
      })))
    }
    rownames(out) <- NULL
    return(out)
  }

  p <- x$indices$overall$p
  if (is.null(p) || length(p) == 0)
    return(data.frame())
  out <- do.call(rbind, lapply(seq_along(p), function(i) {
    mat <- p[[i]]
    expand <- expand.grid(from_group = rownames(mat), to_group = colnames(mat),
                          stringsAsFactors = FALSE)
    data.frame(band = x$bands[i], expand, exposure = as.vector(mat))
  }))
  rownames(out) <- NULL
  out
}

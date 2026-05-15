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

spseg_dots <- function(...) {
  dots <- list(...)
  if ("tol" %in% names(dots)) {
    message("'tol' is ignored by spseg(); zero entropy terms are handled directly")
    dots$tol <- NULL
  }

  dots
}

spseg_output <- function(output) {
  match.arg(output, c("legacy", "full", "indices", "localenv"))
}

spseg_comparison <- function(comparison) {
  match.arg(comparison, c("overall", "pairwise", "both"))
}

spseg_comparison_flags <- function(comparison) {
  comparison <- spseg_comparison(comparison)
  c(overall = comparison %in% c("overall", "both"),
    pairwise = comparison %in% c("pairwise", "both"))
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
                         output, call) {
  keep_inputs <- !identical(output, "indices")
  SegResult(
    coords = if (keep_inputs) coords else NULL,
    data = if (keep_inputs) data else NULL,
    env = env,
    bands = bands,
    indices = indices,
    measures = measures,
    weighting = weighting,
    power = power,
    normalize = normalize,
    neighbors = list(type = "radius", values = bands,
                     comparison = comparison),
    proj4string = proj4string,
    call = call
  )
}

spseg_smoothing_config <- function(smoothing, dots) {
  old_names <- c("nrow", "ncol", "nocol", "window", "sigma")
  supplied_old <- intersect(names(dots), old_names)
  deprecated <- character()

  if (is.null(smoothing)) {
    config <- list(smoothing = "none")
  } else if (is.list(smoothing)) {
    config <- smoothing
  } else {
    config <- list(smoothing = smoothing)
    deprecated <- c(deprecated, "smoothing")
  }

  if (length(supplied_old) > 0) {
    deprecated <- c(deprecated, supplied_old)
    for (nm in supplied_old) {
      config[[nm]] <- dots[[nm]]
      dots[[nm]] <- NULL
    }
  }

  if (length(deprecated) > 0)
    warning("Deprecated smoothing argument(s) ignored as top-level inputs: ",
            paste(unique(deprecated), collapse = ", "),
            ". Use smoothing = list(...) instead.", call. = FALSE)

  if (is.null(config$smoothing))
    config$smoothing <- "none"
  if (!is.null(config$nocol) && is.null(config$ncol))
    config$ncol <- config$nocol
  config$nocol <- NULL
  if (is.null(config$nrow))
    config$nrow <- 100
  if (is.null(config$ncol))
    config$ncol <- 100

  list(config = config, dots = dots)
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

spseg_xlogx <- function(x, base) {
  out <- numeric(length(x))
  positive <- x > 0
  out[positive] <- x[positive] * log(x[positive], base = base)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}

spseg_compute_c <- function(data, env, measures) {
  m <- ncol(data)
  measure_flags <- c("exposure" %in% measures, "information" %in% measures,
                     "diversity" %in% measures, "dissimilarity" %in% measures)
  tmp <- .Call("spsegIDX", as.vector(data), as.vector(env), as.integer(m),
               as.integer(measure_flags))
  results <- spseg_empty_results()
  n <- m^2

  if (!is.na(tmp[1])) {
    results$p <- matrix(tmp[1:n], ncol = m, byrow = TRUE)
    rownames(results$p) <- colnames(results$p) <- colnames(data)
  }
  if (!is.na(tmp[n + 1])) results$h <- tmp[n + 1]
  if (!is.na(tmp[n + 2])) results$r <- tmp[n + 2]
  if (!is.na(tmp[n + 3])) results$d <- tmp[n + 3]

  results
}

spseg_compute_r <- function(data, env, measures) {
  results <- spseg_empty_results()
  m <- ncol(data)
  ptsSum <- sum(data)
  ptsRowSum <- rowSums(data)
  ptsColSum <- colSums(data)
  ptsProp <- ptsColSum / ptsSum
  envProp <- env / rowSums(env)

  if ("exposure" %in% measures) {
    P <- matrix(0, nrow = m, ncol = m)
    rownames(P) <- colnames(P) <- colnames(data)
    for (i in 1:m) {
      A <- data[, i] / ptsColSum[i]
      for (j in 1:m)
        P[i, j] <- sum(A * envProp[, j])
    }
    results$p <- P
  }

  if ("information" %in% measures) {
    Ep <- -rowSums(spseg_xlogx(envProp, base = m))
    E <- -sum(spseg_xlogx(ptsProp, base = m))
    results$h <- 1 - (sum(ptsRowSum * Ep) / (ptsSum * E))
  }

  if ("diversity" %in% measures) {
    Ip <- rowSums(envProp * (1 - envProp))
    I <- sum(ptsProp * (1 - ptsProp))
    results$r <- 1 - sum((ptsRowSum * Ip) / (ptsSum * I))
  }

  if ("dissimilarity" %in% measures) {
    I <- sum(ptsProp * (1 - ptsProp))
    constant <- ptsRowSum / (2 * ptsSum * I)
    Dp <- abs(sweep(envProp, 2, ptsProp))
    results$d <- sum(colSums(Dp * constant))
  }

  results
}

spseg_pairwise_compute_r <- function(data, env, measures) {
  measures <- spseg_measures(measures)
  m <- ncol(data)
  empty <- matrix(NA_real_, nrow = m, ncol = m,
                  dimnames = list(colnames(data), colnames(data)))
  out <- list(d = list(), r = list(), h = list())
  if (m < 2)
    return(out)

  dmat <- rmat <- hmat <- empty
  for (a in seq_len(m - 1)) {
    for (b in (a + 1):m) {
      pair_total <- sum(data[, a] + data[, b])
      if (sum(data[, a]) <= 0 || sum(data[, b]) <= 0)
        next
      P <- c(sum(data[, a]), sum(data[, b])) / pair_total
      E <- -sum(spseg_xlogx(P, base = 2))
      I <- sum(P * (1 - P))
      D_acc <- R_acc <- H_acc <- 0

      for (i in seq_len(nrow(data))) {
        T_iab <- data[i, a] + data[i, b]
        L_ab <- env[i, a] + env[i, b]
        if (T_iab <= 0 || L_ab <= 0)
          next
        p <- c(env[i, a], env[i, b]) / L_ab
        if ("dissimilarity" %in% measures)
          D_acc <- D_acc + (T_iab / (2 * pair_total * I)) *
            sum(abs(p - P))
        if ("diversity" %in% measures)
          R_acc <- R_acc + T_iab * sum(p * (1 - p))
        if ("information" %in% measures)
          H_acc <- H_acc - T_iab * sum(spseg_xlogx(p, base = 2))
      }

      if ("dissimilarity" %in% measures)
        dmat[a, b] <- dmat[b, a] <- D_acc
      if ("diversity" %in% measures)
        rmat[a, b] <- rmat[b, a] <- 1 - R_acc / (pair_total * I)
      if ("information" %in% measures)
        hmat[a, b] <- hmat[b, a] <- 1 - H_acc / (pair_total * E)
    }
  }

  if ("dissimilarity" %in% measures)
    out$d <- list(dmat)
  if ("diversity" %in% measures)
    out$r <- list(rmat)
  if ("information" %in% measures)
    out$h <- list(hmat)
  out
}

spseg_from_localenv <- function(env, measures, useC, negative.rm) {
  measures <- spseg_measures(measures)
  prepared <- spseg_prepare_localenv(env, negative.rm)
  results <- if (useC)
    spseg_compute_c(prepared$data, prepared$env, measures)
  else
    spseg_compute_r(prepared$data, prepared$env, measures)

  SegSpatial(results$d, results$r, results$h, results$p,
             env@coords, env@data, env@env, env@proj4string)
}

spseg_indices_from_localenv <- function(env, measures, useC, negative.rm,
                                        comparison = "overall") {
  measures <- spseg_measures(measures)
  comparison_flags <- spseg_comparison_flags(comparison)
  prepared <- spseg_prepare_localenv(env, negative.rm)
  overall <- list(d = numeric(), r = numeric(), h = numeric(), p = list())
  pairwise <- list(d = list(), r = list(), h = list())

  if (comparison_flags["overall"]) {
    overall <- if (useC)
      spseg_compute_c(prepared$data, prepared$env, measures)
    else
      spseg_compute_r(prepared$data, prepared$env, measures)
    if (length(overall$p) > 0)
      overall$p <- list(overall$p)
    else
      overall$p <- list()
  }

  if (comparison_flags["pairwise"])
    pairwise <- spseg_pairwise_compute_r(prepared$data, prepared$env, measures)

  list(overall = overall, pairwise = pairwise)
}

spseg_surface <- function(x, coords, data, smoothing, verbose) {
  nrow <- smoothing$nrow
  ncol <- smoothing$ncol
  window <- smoothing$window
  sigma <- smoothing$sigma
  smoothing <- smoothing$smoothing
  smoothing <- match.arg(smoothing, c("none", "kernel", "equal"),
                         several.ok = FALSE)
  if (smoothing == "none")
    return(list(coords = coords, data = data))

  if (smoothing == "equal")
    return(surface.equal(x, data, nrow, ncol, verbose))

  if (is.null(window)) {
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    window <- matrix(c(x_range[1], y_range[1],
                       x_range[1], y_range[2],
                       x_range[2], y_range[2],
                       x_range[2], y_range[1]),
                     ncol = 2, byrow = TRUE)
  }

  if (is.null(sigma))
    sigma <- min(bw.nrd(coords[, 1]), bw.nrd(coords[, 2]))

  surface.kernel(coords, data, sigma, nrow, ncol, window, verbose)
}

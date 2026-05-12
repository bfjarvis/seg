# ------------------------------------------------------------------------------
# Internal helpers for spseg()
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg_methods <- function(method) {
  method <- match.arg(method, c("exposure", "information", "diversity",
                                "dissimilarity", "all"), several.ok = TRUE)
  if ("all" %in% method)
    method <- c("exposure", "information", "diversity", "dissimilarity")

  method
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

spseg_compute_c <- function(data, env, method) {
  m <- ncol(data)
  method_flags <- c("exposure" %in% method, "information" %in% method,
                    "diversity" %in% method, "dissimilarity" %in% method)
  tmp <- .Call("spsegIDX", as.vector(data), as.vector(env), as.integer(m),
               as.integer(method_flags))
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

spseg_compute_r <- function(data, env, method) {
  results <- spseg_empty_results()
  m <- ncol(data)
  ptsSum <- sum(data)
  ptsRowSum <- rowSums(data)
  ptsColSum <- colSums(data)
  ptsProp <- ptsColSum / ptsSum
  envProp <- env / rowSums(env)

  if ("exposure" %in% method) {
    P <- matrix(0, nrow = m, ncol = m)
    rownames(P) <- colnames(P) <- colnames(data)
    for (i in 1:m) {
      A <- data[, i] / ptsColSum[i]
      for (j in 1:m)
        P[i, j] <- sum(A * envProp[, j])
    }
    results$p <- P
  }

  if ("information" %in% method) {
    Ep <- -rowSums(spseg_xlogx(envProp, base = m))
    E <- -sum(spseg_xlogx(ptsProp, base = m))
    results$h <- 1 - (sum(ptsRowSum * Ep) / (ptsSum * E))
  }

  if ("diversity" %in% method) {
    Ip <- rowSums(envProp * (1 - envProp))
    I <- sum(ptsProp * (1 - ptsProp))
    results$r <- 1 - sum((ptsRowSum * Ip) / (ptsSum * I))
  }

  if ("dissimilarity" %in% method) {
    I <- sum(ptsProp * (1 - ptsProp))
    constant <- ptsRowSum / (2 * ptsSum * I)
    Dp <- abs(envProp - ptsProp)
    results$d <- sum(colSums(Dp * constant))
  }

  results
}

spseg_from_localenv <- function(env, method, useC, negative.rm) {
  method <- spseg_methods(method)
  prepared <- spseg_prepare_localenv(env, negative.rm)
  results <- if (useC)
    spseg_compute_c(prepared$data, prepared$env, method)
  else
    spseg_compute_r(prepared$data, prepared$env, method)

  SegSpatial(results$d, results$r, results$h, results$p,
             env@coords, env@data, env@env, env@proj4string)
}

spseg_surface <- function(x, coords, data, smoothing, nrow, ncol, window,
                          sigma, verbose) {
  smoothing <- match.arg(smoothing, c("none", "kernel", "equal"),
                         several.ok = FALSE)
  if (smoothing == "none")
    return(list(coords = coords, data = data))

  if (smoothing == "equal")
    return(surface.equal(x, data, nrow, ncol, verbose))

  if (missing(window)) {
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    window <- matrix(c(x_range[1], y_range[1],
                       x_range[1], y_range[2],
                       x_range[2], y_range[2],
                       x_range[2], y_range[1]),
                     ncol = 2, byrow = TRUE)
  }

  if (missing(sigma))
    sigma <- min(bw.nrd(coords[, 1]), bw.nrd(coords[, 2]))

  surface.kernel(coords, data, sigma, nrow, ncol, window, verbose)
}

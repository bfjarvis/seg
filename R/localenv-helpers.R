# ------------------------------------------------------------------------------
# Internal helpers for local environment calculations
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
.geometric_mean_distance <- function(x, max_points = 500) {
  if (inherits(x, "dist")) {
    d <- as.numeric(x)
  } else {
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
  }

  d <- d[d > 0]
  if (length(d) == 0)
    return(0)

  exp(mean(log(d)))
}

.restore_env_dimnames <- function(env, data) {
  colnames(env) <- colnames(data)
  rownames(env) <- rownames(data)
  env
}

.localenv_weighting_id <- function(weighting) {
  match(weighting, c("unweighted", "biweight", "inverse", "exponential")) - 1L
}

.localenv_legacy_message <- function(useExp_supplied, scale_supplied) {
  old_args <- c(if (useExp_supplied) "useExp", if (scale_supplied) "scale")
  if (length(old_args) > 0)
    message("Ignoring deprecated localenv argument(s): ",
            paste(old_args, collapse = ", "),
            ". See ?localenv for the current weighting syntax.")
}

.localenv_settings_message <- function(weighting, power, normalize) {
  message("localenv weighting: ", weighting, "; power: ", power,
          "; normalize: ", normalize)
}

.localenv_dots <- function(...) {
  dots <- list(...)
  if ("tol" %in% names(dots)) {
    message("Ignoring deprecated localenv argument: tol. See ?localenv.")
    dots$tol <- NULL
  }
  if (length(dots) > 0)
    stop("unused argument(s): ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  invisible(NULL)
}

localenv_coords <- function(coords, data, power, weighting, normalize, maxdist) {
  if (nrow(coords) != nrow(data))
    stop("'data' must have the same number of rows as 'coords'", call. = FALSE)

  xval <- coords[, 1]
  yval <- coords[, 2]
  dim <- ncol(data)
  env <- .Call("envconstruct", xval, yval, as.vector(data), as.integer(dim),
               power, as.integer(.localenv_weighting_id(weighting)),
               as.integer(normalize), maxdist)

  .restore_env_dimnames(env, data)
}

localenv_dist <- function(sprel, data, power, weighting, normalize, maxdist) {
  sprel <- as.matrix(sprel)
  if (nrow(sprel) != nrow(data))
    stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)

  env <- matrix(nrow = nrow(data), ncol = ncol(data))
  for (i in 1:nrow(data)) {
    in_band <- if (maxdist >= 0) sprel[i, ] <= maxdist else
      rep(TRUE, nrow(sprel))
    bw <- if (normalize) maxdist else 1
    if (normalize && maxdist < 0)
      bw <- .geometric_mean_distance(as.dist(sprel))
    if (bw <= 0)
      bw <- 1

    d <- sprel[i, ]
    positive_neighbors <- d[in_band & d > 0]
    local_tol <- if (length(positive_neighbors) > 0) min(positive_neighbors) / 2 else
      bw
    d_inverse <- d
    d_inverse[d_inverse == 0] <- local_tol

    weight <- switch(weighting,
      unweighted = rep(1, length(d)),
      biweight = (1 - (d / bw)^power)^power,
      inverse = 1 / ((d_inverse / bw)^power),
      exponential = exp(-d / bw)
    )

    if (weighting == "biweight")
      weight[d > bw] <- 0
    weight[!in_band] <- 0
    weight[weight < 0] <- 0
    env[i, ] <- apply(data, 2, function(z) sum(z * weight) / sum(weight))
  }

  .restore_env_dimnames(env, data)
}

localenv_nb <- function(sprel, data) {
  xmat <- spdep::nb2mat(sprel, style = "W")
  if (nrow(xmat) != nrow(data))
    stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)

  env <- xmat %*% data
  .restore_env_dimnames(env, data)
}

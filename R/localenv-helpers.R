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

localenv_coords <- function(coords, data, power, useExp, scale, maxdist, tol) {
  if (nrow(coords) != nrow(data))
    stop("'data' must have the same number of rows as 'coords'", call. = FALSE)

  xval <- coords[, 1]
  yval <- coords[, 2]
  dim <- ncol(data)
  env <- .Call("envconstruct", xval, yval, as.vector(data), as.integer(dim),
               power, as.integer(useExp), as.integer(scale), maxdist, tol)

  .restore_env_dimnames(env, data)
}

localenv_dist <- function(sprel, data, power, useExp, scale, maxdist, tol) {
  sprel <- as.matrix(sprel)
  if (nrow(sprel) != nrow(data))
    stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)

  env <- matrix(nrow = nrow(data), ncol = ncol(data))
  for (i in 1:nrow(data)) {
    if (scale) {
      d <- sprel[i, ] / maxdist
      if (power == 0)
        weight <- rep(1, length(d))
      else if (useExp)
        weight <- (exp(d * power * -1) - exp(power * -1)) /
          (1 - exp(power * -1))
      else
        weight <- (1 - d^power)^power
    } else if (useExp) {
      weight <- exp(power * sprel[i, ] * -1)
    } else {
      weight <- 1 / (sprel[i, ] + tol)^power
    }

    if (maxdist >= 0)
      weight[which(sprel[i, ] > maxdist)] <- 0
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

# ------------------------------------------------------------------------------
# Function 'localenv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
.geometric_mean_distance <- function(x, max_points = 500) {
  if (inherits(x, "dist"))
    d <- as.numeric(x)
  else {
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

localenv <- function(x, data, power = 2, useExp = TRUE, scale = TRUE,
                     maxdist, sprel, tol = .Machine$double.eps) {
  
  tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; data <- tmp$data; proj4string <- tmp$proj4string
  
  maxdist_missing <- missing(maxdist)
  if (maxdist_missing)
    maxdist <- -1
  else if (!is.numeric(maxdist))
    stop("'maxdist' must be numeric", call. = FALSE)
  else if (maxdist < 0)
    stop("'maxdist' must be greater than or equal to 0", call. = FALSE)
  
  if (missing(sprel)) 
    sprel <- coords
  else if (!inherits(sprel, "nb") && !inherits(sprel, "dist"))
    stop("invalid object 'sprel'", call. = FALSE)

  if (scale && !inherits(sprel, "nb")) {
    if (maxdist_missing) {
      maxdist <- .geometric_mean_distance(sprel)
    }
    if (maxdist <= 0)
      stop("'maxdist' must be greater than 0 when 'scale' is TRUE",
           call. = FALSE)
  }
  
  env <- localenv.get(sprel, data, power, useExp, scale, maxdist, tol)
  
  SegLocal(coords, data, env, st_crs(proj4string))
}

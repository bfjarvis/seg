# ------------------------------------------------------------------------------
# Function 'localenv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
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
  
  if (missing(sprel)) {
    relation <- coords
    relation_type <- "coords"
  } else if (inherits(sprel, "dist")) {
    relation <- sprel
    relation_type <- "dist"
  } else if (inherits(sprel, "nb")) {
    relation <- sprel
    relation_type <- "nb"
  } else {
    stop("invalid object 'sprel'", call. = FALSE)
  }

  if (scale && relation_type != "nb") {
    if (maxdist_missing) {
      maxdist <- .geometric_mean_distance(relation)
    }
    if (maxdist <= 0)
      stop("'maxdist' must be greater than 0 when 'scale' is TRUE",
           call. = FALSE)
  }

  env <- switch(relation_type,
    coords = localenv_coords(relation, data, power, useExp, scale, maxdist, tol),
    dist = localenv_dist(relation, data, power, useExp, scale, maxdist, tol),
    nb = localenv_nb(relation, data)
  )
  
  SegLocal(coords, data, env, st_crs(proj4string))
}

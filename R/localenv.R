# ------------------------------------------------------------------------------
# Function 'localenv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
localenv <- function(x, data, power = 3,
                     weighting = c("biweight", "unweighted", "inverse",
                                   "exponential"),
                     normalize = TRUE, maxdist, sprel,
                     useExp = NULL, scale = NULL, ...) {
  
  tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; data <- tmp$data; proj4string <- tmp$proj4string
  .localenv_dots(...)
  useExp_supplied <- !missing(useExp)
  scale_supplied <- !missing(scale)
  if (!useExp_supplied)
    useExp <- NULL
  if (!scale_supplied)
    scale <- NULL
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)
  .localenv_legacy_message(useExp_supplied, scale_supplied)
  .localenv_settings_message(weighting, power, normalize)
  
  maxdist_missing <- missing(maxdist)
  if (maxdist_missing)
    maxdist <- -1
  else if (!is.numeric(maxdist))
    stop("'maxdist' must be numeric", call. = FALSE)
  else if (maxdist < 0)
    stop("'maxdist' must be greater than or equal to 0", call. = FALSE)

  if (!maxdist_missing && maxdist == 0) {
    message("localenv maxdist is 0; using observed data as local environments")
    return(SegLocal(coords, data, data, st_crs(proj4string)))
  }
  
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

  if (normalize && relation_type != "nb") {
    if (maxdist_missing) {
      maxdist <- .geometric_mean_distance(relation)
    }
    if (maxdist <= 0)
      stop("'maxdist' must be greater than 0 when 'normalize' is TRUE",
           call. = FALSE)
  }

  env <- switch(relation_type,
    coords = localenv_coords(relation, data, power, weighting, normalize,
                             maxdist),
    dist = localenv_dist(relation, data, power, weighting, normalize,
                         maxdist),
    nb = localenv_nb(relation, data)
  )
  
  SegLocal(coords, data, env, st_crs(proj4string))
}

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

  dots <- list(...)
  if ("tol" %in% names(dots)) {
    message("Ignoring deprecated localenv argument: tol. See ?localenv.")
    dots$tol <- NULL
  }
  if (length(dots) > 0)
    stop("unused argument(s): ", paste(names(dots), collapse = ", "),
         call. = FALSE)

  if (!missing(sprel))
    stop("'sprel' is no longer supported. Distance and neighbor-list local ",
         "environment inputs have been deprecated; use coordinates or sf ",
         "geometry instead.", call. = FALSE)

  useExp_supplied <- !missing(useExp)
  scale_supplied <- !missing(scale)
  if (!useExp_supplied)
    useExp <- NULL
  if (!scale_supplied)
    scale <- NULL
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)
  old_args <- c(if (useExp_supplied) "useExp", if (scale_supplied) "scale")
  if (length(old_args) > 0)
    message("Ignoring deprecated localenv argument(s): ",
            paste(old_args, collapse = ", "),
            ". See ?localenv for the current weighting syntax.")
  message("localenv weighting: ", weighting, "; power: ", power,
          "; normalize: ", normalize)

  maxdist_missing <- missing(maxdist)
  if (!maxdist_missing) {
    if (!is.numeric(maxdist))
      stop("'maxdist' must be numeric", call. = FALSE)
    if (maxdist < 0)
      stop("'maxdist' must be greater than or equal to 0", call. = FALSE)
  }

  if (maxdist_missing) {
    tmp <- if (missing(data))
      suppressMessages(chksegdata(x))
    else
      suppressMessages(chksegdata(x, data))
    maxdist <- if (normalize) {
      .geometric_mean_distance(tmp$coords)
    } else {
      d <- as.numeric(dist(tmp$coords))
      if (length(d) > 0) max(d) else 0
    }
    if (normalize && maxdist <= 0)
      stop("'maxdist' must be greater than 0 when 'normalize' is TRUE",
           call. = FALSE)
  }

  if (!maxdist_missing && maxdist == 0)
    message("localenv maxdist is 0; using observed data as local environments")

  args <- list(
    x = x,
    maxdist = maxdist,
    power = power,
    weighting = weighting,
    normalize = normalize,
    output = "localenv"
  )
  if (!missing(data))
    args$data <- data

  result <- do.call(spseg, args)
  SegLocal(result@coords, result@data, result@env[[1]], result@proj4string,
           geometry = result@geometry)
}

# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- function(x,
                  data,
                  measures = "all",
                  sprel = NULL,
                  maxdist = NULL,
                  bands = NULL,
                  weighting = c("biweight", "unweighted", "inverse", "exponential"),
                  power = 3,
                  normalize = TRUE,
                  smoothing = NULL,
                  comparison = c("overall", "pairwise", "both"),
                  output = c("legacy", "full", "indices", "localenv"),
                  useC = TRUE,
                  negative.rm = FALSE,
                  verbose = FALSE,
                  ...) {

  call <- match.call()
  output <- spseg_output(output)
  comparison <- spseg_comparison(comparison)

  # parse inputs
  dots <- spseg_dots(...)
  if ("method" %in% names(dots)) {
    warning("'method' is deprecated; use 'measures' instead.",
            call. = FALSE)
    if (missing(measures))
      measures <- dots$method
    dots$method <- NULL
  }
  smoothing_config <- spseg_smoothing_config(smoothing, dots)
  dots <- smoothing_config$dots
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)
  if (!identical(comparison, "overall") && identical(output, "legacy"))
    output <- "full"

  data_arg <- if (missing(data)) NULL else data
  bands <- spseg_bands(maxdist, bands, output, x, data_arg)
  if (!is.null(bands) && any(!is.finite(bands)))
    stop("'bands' must contain finite numeric values", call. = FALSE)
  if (!is.null(bands) && any(bands < 0))
    stop("'bands' must be greater than or equal to 0", call. = FALSE)
  if (!is.null(bands) && length(bands) > 1 && identical(output, "legacy"))
    output <- "full"

  localenv_args <- list(power = power, weighting = weighting,
                        normalize = normalize)
  if (!is.null(bands) && length(bands) == 1)
    localenv_args$maxdist <- bands
  if (!is.null(sprel))
    localenv_args$sprel <- sprel

  # Existing one-band workflows keep their legacy return by default.
  if (inherits(x, "SegLocal") && identical(output, "legacy"))
    return(spseg_from_localenv(x, measures, useC, negative.rm))

  if (inherits(x, "SegLocal")) {
    indices <- if (identical(output, "localenv")) list() else
      spseg_indices_from_localenv(x, measures, useC, negative.rm,
                                  comparison = comparison)
    if (length(indices) > 0)
      indices <- spseg_indices_from_engine(indices, 1)
    env <- if (identical(output, "indices")) NULL else list(x@env)
    return(spseg_result(
      coords = x@coords, data = x@data, env = env, bands = 1,
      indices = indices,
      measures = if (length(indices) > 0) spseg_measures(measures) else
        character(),
      comparison = comparison, weighting = character(), power = numeric(),
      normalize = logical(),
      proj4string = x@proj4string, output = output, call = call,
      geometry = x@geometry
    ))
  }

  # verify data and prepare the analysis surface
  checked <- if (verbose) chksegdata(x, data) else
    suppressMessages(chksegdata(x, data))
  surface <- if (identical(smoothing_config$config$smoothing, "none")) {
    list(coords = checked$coords, data = checked$data)
  } else {
    spseg_surface(
      x,
      checked$coords,
      checked$data,
      smoothing_config$config,
      verbose
    )
  }

  if (!is.null(bands) && length(bands) > 1 && !is.null(sprel))
    stop("multiple bands are not yet supported with 'sprel'", call. = FALSE)

  if (!is.null(sprel) && !identical(output, "legacy")) {
    env <- do.call(
      localenv,
      c(list(x = surface$coords, data = surface$data), localenv_args, dots)
    )
    env <- update(env, proj4string = st_crs(checked$proj4string))
    indices <- if (identical(output, "localenv")) list() else
      spseg_indices_from_localenv(env, measures, useC, negative.rm,
                                  comparison = comparison)
    if (length(indices) > 0)
      indices <- spseg_indices_from_engine(indices, bands)
    return(spseg_result(
      coords = surface$coords, data = surface$data,
      env = if (identical(output, "indices")) NULL else list(env@env),
      bands = bands, indices = indices,
      measures = if (length(indices) > 0) spseg_measures(measures) else
        character(),
      comparison = comparison, weighting = weighting, power = power,
      normalize = normalize,
      proj4string = st_crs(checked$proj4string), output = output,
      call = call, geometry = checked$geometry
    ))
  }

  # New unified result path: multiple bands, indices-only, localenv-only, or full.
  if (!is.null(bands) && (length(bands) > 1 || !identical(output, "legacy"))) {
    measures_to_compute <- if (identical(output, "localenv")) character() else
      measures
    engine <- seg_engine_coords(
      coords = surface$coords,
      data = surface$data,
      bands = bands,
      power = power,
      weighting = weighting,
      normalize = normalize,
      measures = measures_to_compute,
      comparison = comparison,
      keep_env = !identical(output, "indices"),
      keep_indices = !identical(output, "localenv")
    )
    indices <- spseg_indices_from_engine(engine$indices, bands)
    return(spseg_result(
      coords = surface$coords, data = surface$data, env = engine$env,
      bands = bands, indices = indices,
      measures = if (length(indices) > 0) spseg_measures(measures) else
        character(),
      comparison = comparison, weighting = weighting, power = power,
      normalize = normalize,
      proj4string = st_crs(checked$proj4string), output = output,
      call = call, geometry = checked$geometry
    ))
  }

  # Legacy path: construct one local environment object, then calculate measures.
  env <- do.call(
    localenv,
    c(list(x = surface$coords, data = surface$data), localenv_args, dots)
  )
  env <- update(env, proj4string = st_crs(checked$proj4string))

  spseg_from_localenv(env, measures, useC, negative.rm)
}

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
                  neighbors = c("radius", "knn"),
                  search = c("kdtree", "brute"),
                  surface = c("raw", "grid", "pycno"),
                  smoothing = NULL,
                  comparison = c("overall", "pairwise", "both"),
                  output = c("indices", "full", "localenv", "legacy"),
                  useC = TRUE,
                  negative.rm = FALSE,
                  verbose = FALSE,
                  ...) {

  call <- match.call()
  output <- spseg_output(output)
  comparison <- spseg_comparison(comparison)
  neighbors <- spseg_neighbors(neighbors)
  search <- spseg_search(search)
  surface <- spseg_surface_method(surface)

  # parse inputs
  dots <- spseg_dots(...)
  if ("method" %in% names(dots)) {
    warning("'method' is deprecated; use 'measures' instead.",
            call. = FALSE)
    if (missing(measures))
      measures <- dots$method
    dots$method <- NULL
  }
  surface_config <- spseg_surface_config(surface, smoothing, dots)
  surface <- surface_config$surface
  dots <- surface_config$dots
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)
  data_arg <- if (missing(data)) NULL else data
  if (!is.null(sprel))
    stop("'sprel' is no longer supported. Distance and neighbor-list local ",
         "environment inputs have been deprecated; use coordinates or sf ",
         "geometry instead.", call. = FALSE)
  if (neighbors == "knn" && !is.null(maxdist))
    stop("'maxdist' is not supported with kNN neighborhoods; use 'bands' for population thresholds",
         call. = FALSE)
  if (neighbors == "knn" && is.null(bands))
    stop("'bands' must be supplied for count-based kNN neighborhoods",
         call. = FALSE)

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
      geometry = x@geometry, neighbors = "radius", search = "brute"
    ))
  }

  bands <- spseg_bands(maxdist, bands, output, x, data_arg)
  if (!is.null(bands) && any(!is.finite(bands)))
    stop("'bands' must contain finite numeric values", call. = FALSE)
  if (!is.null(bands) && neighbors == "radius" && any(bands < 0))
    stop("'bands' must be greater than or equal to 0", call. = FALSE)
  if (!is.null(bands) && neighbors == "knn" && any(bands <= 0))
    stop("'bands' must be greater than 0 for kNN neighborhoods",
         call. = FALSE)
  if (!is.null(bands) && length(bands) > 1 && identical(output, "legacy"))
    stop("'output = \"legacy\"' supports only one bandwidth", call. = FALSE)

  # verify data and prepare the analysis surface
  checked <- if (verbose) chksegdata(x, data) else
    suppressMessages(chksegdata(x, data))
  pop_surface <- spseg_surface(
    x = x,
    coords = checked$coords,
    data = checked$data,
    surface = surface,
    args = dots,
    verbose = verbose
  )
  pop_surface$coords <- as.matrix(pop_surface$coords)
  pop_surface$data <- as.matrix(pop_surface$data)

  if (is.null(bands)) {
    bands <- if (identical(output, "legacy")) {
      if (normalize) .geometric_mean_distance(pop_surface$coords) else -1
    } else {
      default_bands(pop_surface$coords, pop_surface$data)
    }
    if (normalize && bands <= 0)
      stop("'maxdist' must be greater than 0 when 'normalize' is TRUE",
           call. = FALSE)
  }

  measures_to_compute <- if (identical(output, "localenv")) character() else
    measures
  engine <- seg_engine_coords(
    coords = pop_surface$coords,
    data = pop_surface$data,
    bands = bands,
    power = power,
    weighting = weighting,
    normalize = normalize,
    measures = measures_to_compute,
    comparison = comparison,
    neighbors = neighbors,
    search = search,
    keep_env = !identical(output, "indices"),
    keep_indices = !identical(output, "localenv")
  )
  indices <- spseg_indices_from_engine(engine$indices, bands)
  result <- spseg_result(
    coords = pop_surface$coords, data = pop_surface$data, env = engine$env,
    bands = bands, indices = indices,
    measures = if (!identical(output, "localenv")) spseg_measures(measures) else
      character(),
    comparison = comparison, weighting = weighting, power = power,
    normalize = normalize,
    proj4string = st_crs(checked$proj4string), output = output,
    call = call,
    geometry = if (identical(surface, "raw")) checked$geometry else
      pop_surface$geometry,
    neighbors = neighbors,
    search = search
  )
  if (identical(output, "legacy"))
    return(spseg_legacy_from_result(result))
  result
}

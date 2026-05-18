# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- function(x,
                  data = NULL,
                  measures = "all",
                  bands = NULL,
                  weighting = c("biweight", "unweighted", "inverse", "exponential"),
                  power = 3,
                  normalize = TRUE,
                  neighbors = c("radius", "knn"),
                  search = c("kdtree", "brute"),
                  surface = c("raw", "grid", "pycno"),
                  comparison = c("overall", "pairwise", "both"),
                  output = c("indices", "full", "localenv"),
                  verbose = FALSE,
                  ...) {

  call <- match.call()
  output <- spseg_output(output)
  comparison <- spseg_comparison(comparison)
  neighbors <- spseg_neighbors(neighbors)
  search <- spseg_search(search)
  surface <- spseg_surface_method(surface)
  dots <- list(...)
  weighting <- if (missing(weighting)) "biweight" else match.arg(weighting)
  normalize <- isTRUE(normalize)

  if (neighbors == "knn" && is.null(bands))
    stop("'bands' must be supplied for count-based kNN neighborhoods",
         call. = FALSE)

  bands <- spseg_bands(bands)
  if (!is.null(bands) && any(!is.finite(bands)))
    stop("'bands' must contain finite numeric values", call. = FALSE)
  if (!is.null(bands) && neighbors == "radius" && any(bands < 0))
    stop("'bands' must be greater than or equal to 0", call. = FALSE)
  if (!is.null(bands) && neighbors == "knn" && any(bands <= 0))
    stop("'bands' must be greater than 0 for kNN neighborhoods",
         call. = FALSE)

  checked <- if (verbose) {
    if (is.null(data)) chksegdata(x) else chksegdata(x, data)
  } else {
    if (is.null(data)) suppressMessages(chksegdata(x)) else
      suppressMessages(chksegdata(x, data))
  }
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
    bands <- default_bands(pop_surface$coords, pop_surface$data)
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
  indices <- if (identical(output, "localenv")) {
    list()
  } else {
    spseg_indices_from_engine(engine$indices, bands)
  }
  result <- spseg_result(
    coords = pop_surface$coords, data = pop_surface$data, env = engine$env,
    bands = bands, indices = indices,
    measures = if (!identical(output, "localenv")) spseg_measures(measures) else
      character(),
    comparison = comparison, weighting = weighting, power = power,
    normalize = normalize,
    crs = st_crs(checked$proj4string), output = output,
    call = call,
    geometry = if (identical(surface, "raw")) checked$geometry else
      pop_surface$geometry,
    neighbors = neighbors,
    search = search,
    surface = surface
  )
  result
}

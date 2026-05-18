#' Index of Spatial Proximity
#'
#' Computes the index of spatial proximity developed by White (1983). This
#' measure estimates the level of clustering by comparing the average distance
#' between members of the same group with that between all individuals
#' (regardless of the groups to which they belong). The results may change
#' drastically depending on the definition of distance.
#'
#' `nb` must be a square matrix (i.e., one that has the same number of rows
#' and columns) but does not have to be symmetric. If `nb` is not given,
#' `isp` attempts to create a distance matrix of `x` using the function
#' [stats::dist()] and use it as `nb`. The optional arguments in `...` will
#' be passed to [stats::dist()].
#'
#' @param x A numeric matrix or data frame with coordinates (each row is a
#'   point), or an object of class [sf::sf].
#' @param data An object of class `matrix`, or one that can be coerced to
#'   that class. The number of rows in `data` should equal the number of
#'   geographic units in `x`, and the number of columns should be greater
#'   than one (i.e., at least two population groups are required). This can
#'   be missing if `x` has a data frame attached to it.
#' @param nb An optional `matrix` object indicating the distances between
#'   the geographic units.
#' @param fun A function for the calculation of proximity. The function
#'   should take a numeric vector as an argument (distance) and return a
#'   vector of the same length (proximity). If this is not specified, a
#'   negative exponential function is used by default.
#' @param verbose Logical. If TRUE, print the current stage of the
#'   computation and time spent on each job to the screen.
#' @param ... Optional arguments to be passed to [stats::dist()] when
#'   calculating the distances between the geographic units in `x`. Ignored
#'   if `nb` is given. See [stats::dist()] for available options.
#'
#' @return A single numeric value indicating the degree of segregation; a
#'   value of 1 indicates absence of segregation, and values greater than 1.0
#'   indicate clustering. If the index value is less than one, it indicates
#'   an unusual form of segregation (i.e., people live closer to other
#'   population groups).
#'
#' @references
#' White, M. J. (1983). The measurement of spatial segregation.
#'   *The American Journal of Sociology*, **88**, 1008-1018.
#'
#' @seealso [dissim()], [stats::dist()]
#'
#' @examples
#' # uses the idealised landscapes in 'segdata'
#' data(segdata)
#' geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,0)))) |>
#'   sf::st_make_grid(cellsize = 1)
#' grd_sf <- sf::st_sf(geom)
#'
#' d <- rep(NA, 8) # index of dissimilarity
#' p <- rep(NA, 8) # index of spatial proximity
#' for (i in 1:8) {
#'   idx <- 2 * i
#'   d[i] <- dissim(data = segdata[,(idx-1):idx])$d
#'   p[i] <- isp(grd_sf, data = segdata[,(idx-1):idx])
#' }
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_data <- do.call(
#'     rbind,
#'     lapply(1:8, function(i) {
#'       idx <- 2 * i - 1
#'       data.frame(
#'         pattern = paste0(
#'           LETTERS[i], ": D = ", round(d[i], 2), ", P = ", round(p[i], 2)
#'         ),
#'         group_share = segdata[, idx],
#'         geometry = sf::st_geometry(grd_sf)
#'       )
#'     })
#'   )
#'   plot_data <- sf::st_sf(plot_data)
#'
#'   ggplot2::ggplot(plot_data) +
#'     ggplot2::geom_sf(ggplot2::aes(fill = group_share), color = "white") +
#'     ggplot2::facet_wrap(~ pattern) +
#'     ggplot2::scale_fill_viridis_c() +
#'     ggplot2::theme_void()
#' }
#'
#' @export
isp <- function(x, data, nb, fun, verbose = FALSE, ...) {
  # The input object 'x' can be either a class of 'sf' or 'data.frame'.
  # Depending on the class of 'x', 'data' may not be required. The internal
  # function 'chksegdata()' processes the information as required by the
  # current function.

  # Process input data using chksegdata()
  if (verbose) {
    tmp <- chksegdata(x, data)
  } else {
    tmp <- suppressMessages(chksegdata(x, data))
  }

  coords <- tmp$coords
  pdf <- tmp$data

  # Verify 'coords' and 'data'
  if (ncol(pdf) < 2) {
    stop("'data' must be a matrix with at least two columns", call. = FALSE)
  } else if (!is.numeric(pdf)) {
    stop("'data' must be a numeric matrix", call. = FALSE)
  } else if (nrow(pdf) != nrow(coords)) {
    stop("'data' must have the same number of rows as 'x'", call. = FALSE)
  }

  # Generate distance matrix if 'nb' is missing
  if (missing(nb)) {
    nb <- as.matrix(dist(coords, ...))
  }

  # Use default negative exponential function if 'fun' is missing
  if (missing(fun)) {
    fun <- function(z) exp(-z)
  }

  # Process distance matrix
  if (isSymmetric(nb)) {
    pairID <- t(combn(1:nrow(nb), 2)) # Unique pairs
    pairDist <- as.numeric(as.dist(nb))
  } else {
    pairID <- expand.grid(1:nrow(nb), 1:nrow(nb))
    pairDist <- as.numeric(nb)
  }

  VALID <- which(pairDist != 0)
  pairID <- pairID[VALID, ]
  pairDist <- pairDist[VALID]

  # Compute spatial interaction effect
  speffect <- fun(pairDist)

  pRow <- rowSums(pdf) # Total population by census tracts
  pCol <- colSums(pdf) # Total population by groups

  pA <- sapply(1:ncol(pdf), function(i) {
    sum(pdf[pairID[, 1], i] * pdf[pairID[, 2], i] * speffect) / pCol[i]
  })

  sum(pA) / (sum(pRow[pairID[, 1]] * pRow[pairID[, 2]] * speffect) / sum(pCol))
}

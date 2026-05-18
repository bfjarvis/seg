#' Index of Dissimilarity
#'
#' Calculates the index of dissimilarity proposed by Duncan and Duncan (1955).
#' If `x` is given and `adjust = TRUE`, the index is adjusted to reflect the
#' spatial distribution of population.
#'
#' `dissim` calculates the index of dissimilarity for `data`. If `data` is
#' missing, it attempts to extract the data from `x`. If `x` is not given, or
#' if it is not a valid [sf::sf] object, the function stops with an error.
#'
#' When `x` is given and `adjust` is set to TRUE, the index is adjusted
#' according to the suggestions of Morrill (1991) and Wong (1993). This
#' automatic adjustment uses [spdep] for Morrill's D(adj) and [terra] for
#' Wong's D(w) and D(s), when those packages are available.
#'
#' @param x An optional [sf::sf] object with polygon geometry.
#' @param data A numeric matrix or data frame with two columns that represent
#'   mutually exclusive population groups (e.g., Asians and non-Asians). If
#'   more than two columns are given, only the first two will be used for
#'   computing the index.
#' @param adjust Logical. If TRUE, and if `x` is given, the index of
#'   dissimilarity is adjusted according to the suggestions of Morrill (1991)
#'   and Wong (1993), depending on packages installed on the system. See
#'   Details for more information. Ignored if `x` is not given.
#' @param use_queen Logical. If TRUE, the queen contiguity will be used.
#' @param verbose Logical. If TRUE, print the current stage of the computation
#'   and time spent on each job to the screen.
#'
#' @return A list containing:
#'   \item{d}{Index of dissimilarity.}
#'   \item{dm}{Index of dissimilarity adjusted according to Morrill (1991).
#'     NA if not calculated.}
#'   \item{dw}{Index of dissimilarity adjusted according to Wong (1993).
#'     NA if not calculated.}
#'   \item{ds}{Index of dissimilarity adjusted according to Wong (1993).
#'     NA if not calculated.}
#'
#' @references
#' Duncan, O. D., & Duncan, B. (1955). A methodological analysis of segregation
#'   indexes. *American Sociological Review*, **20**, 210-217.
#'
#' Morrill, R. L. (1991). On the measure of geographic segregation.
#'   *Geography Research Forum*, **11**, 25-36.
#'
#' Wong, D. W. S. (1993). Spatial indices of segregation. *Urban Studies*,
#'   **30**, 559-572.
#'
#' @seealso [spseg()]
#'
#' @examples
#' if (requireNamespace("spdep", quietly = TRUE)) {
#'   # uses the idealised landscapes in 'segdata'
#'   data(segdata)
#'   geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,0)))) |>
#'     sf::st_make_grid(cellsize = 1)
#'   grd_sf <- sf::st_sf(geom)
#'   grd.nb <- spdep::nb2mat(spdep::poly2nb(grd_sf, queen = FALSE), style = "B")
#'   grd.nb <- grd.nb / sum(grd.nb)
#'
#'   d <- list(); m <- list()
#'   for (i in 1:8) {
#'     idx <- 2 * i
#'     d[[i]] <- dissim(data = segdata[,(idx-1):idx])
#'     m[[i]] <- dissim(grd_sf, segdata[,(idx-1):idx], adjust = TRUE)
#'   }
#' }
#'
#' @export
dissim <- function(
  x,
  data,
  adjust = FALSE,
  use_queen = FALSE,
  verbose = FALSE
) {
  # If a sf object 'x' is provided:
  if (!missing(x)) {
    if (verbose) {
      tmp <- chksegdata(x, data)
    } else {
      tmp <- suppressMessages(chksegdata(x, data))
    }

    coords <- tmp$coords
    data <- tmp$data
  }

  if (ncol(data) > 2) {
    warning(
      "'data' has more than two columns; only the first two are used",
      call. = FALSE
    )
    data <- data[, 1:2]
  }

  # 'data' is supposed to contain population. It cannot be negative values.
  # Makes an error if any of the values is less than 0.
  if (any(data < 0)) {
    stop("negative value(s) in 'data'", call. = FALSE)
  }

  colsum <- apply(data, 2, sum)
  if (any(colsum <= 0)) {
    stop("the sum of each column in 'data' must be > 0", call. = FALSE)
  }

  out <- list(d = NA, dm = NA, dw = NA, ds = NA)

  # ****************************************************************************
  #
  # Duncan and Duncan's index of dissimilarity
  #
  # ****************************************************************************
  b <- data[, 1] / sum(data[, 1]) # Blacks
  w <- data[, 2] / sum(data[, 2]) # Whites
  out$d <- as.vector(sum(abs(b - w)) / 2)

  # ****************************************************************************
  #
  # Adjust the original D index if 'adjust' is TRUE and 'x' is given.
  #
  # The function attempts the following adjustments.
  #   (1) Morrill (1991) if 'spdep' is available,
  #   (2) Wong (1993) if 'spdep' and 'terra' are available,
  #   (3) Custom adjustment if an user defined 'nb' is given.
  #
  # ****************************************************************************
  if (!missing(x) & adjust) {
    userpkg <- .packages(all.available = TRUE) # Find out packages available

    if ("spdep" %in% userpkg) {
      tmp <- tryCatch(
        .use_contiguity(x, data, use_queen, verbose),
        error = function(e) print(e)
      )
      if (is.numeric(tmp)) {
        out$dm <- out$d - tmp
      } else if (verbose) {
        message("failed to calculate D(adj)")
      }

      if ("terra" %in% userpkg) {
        tmp <- tryCatch(
          .use_common_boundary(x, data, verbose),
          error = function(e) print(e)
        )
        if (is.numeric(tmp[1])) {
          out$dw <- out$d - tmp[1]
        } else if (verbose) {
          message("failed to calculate D(w)")
        }
        if (is.numeric(tmp[2])) {
          out$ds <- out$d - tmp[2]
        } else if (verbose) {
          message("failed to calculate D(s)")
        }
      }
    }
  }

  return(out)
}

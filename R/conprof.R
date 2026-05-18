#' Concentration Profile
#'
#' Draws a graph that shows the pattern of residential concentration for a
#' population group and calculates its summary statistic as suggested by
#' Hong and Sadahiro (2013).
#'
#' For `n` equally-spaced thresholds between 0 and 1, `conprof` identifies
#' the areas where the selected population group comprises at least the given
#' threshold proportions; computes how many of the group members live in these
#' areas; and plots them on a 2D plane with the threshold values in the
#' horizontal axis and the proportions of the people in the vertical axis.
#'
#' The summary statistic is calculated by estimating the area between the
#' concentration profile and a hypothetical line that represents a uniform
#' distribution (see the examples).
#'
#' @param data An object of class `matrix`, or one that can be coerced to
#'   that class. Each column represents a population group. The number of
#'   columns should be greater than one (i.e., at least two population
#'   groups are required).
#' @param grpID A numeric value specifying the population group (i.e., column
#'   in `data`) to be analysed. Multiple values are not allowed.
#' @param n A numeric value indicating the number of thresholds to be used.
#'   A large value of `n` creates a smoother-looking graph but slows down
#'   the calculation speed.
#' @param graph Logical. If TRUE, draw the concentration profile for the
#'   specified population group.
#' @param add Logical. If TRUE, add the graph to the current plot.
#' @param ... Optional arguments to be passed to [graphics::plot()] when
#'   `add` is FALSE, or to [graphics::lines()] otherwise. Ignored when
#'   `graph` is FALSE.
#'
#' @return A `list` object with the following three elements:
#'   \item{x}{the threshold values.}
#'   \item{y}{the proportions of the people who live in the areas where they
#'     comprise at least the corresponding threshold percentages in the local
#'     population composition.}
#'   \item{d}{the summary statistic for the concentration profile.}
#'
#' @references
#' Poulsen, M., Johnston, R., and Forrest J. (2002) Plural cities and ethnic
#'   enclaves: Introducing a measurement procedure for comparative study.
#'   *International Journal of Urban and Regional Research*, **26**, 229-243.
#'
#' Hong, S.-Y. and Sadahiro, Y. (2013) Measuring geographic segregation:
#'   A graph-based approach. *Journal of Geographical Systems*, **16**,
#'   211-231.
#'
#' @examples
#' xx <- runif(100) # random distribution
#' xx <- xx * (4000 / sum(xx))
#' yy <- rep(c(40, 60), 100) # no segregation
#' zz <- rep(c(100, 0), c(40, 60)) # complete segregation
#'
#' set1 <- cbind(xx, 100 - xx)
#' set2 <- matrix(yy, ncol = 2, byrow = TRUE)
#' set3 <- cbind(zz, 100 - zz)
#'
#' par(mar = c(5.1, 4.1, 2.1, 2.1))
#' out1 <- conprof(set1, grpID = 1,
#'   xlab = "Threshold level (%)",
#'   ylab = "Population proportion (%)",
#'   cex.lab = 0.9, cex.axis = 0.9, lty = "dotted")
#' out2 <- conprof(set2, grpID = 1, add = TRUE,
#'   lty = "longdash")
#' out3 <- conprof(set3, grpID = 1, add = TRUE)
#' title(main = paste("R =", round(out1$d, 2)))
#'
#' # shaded areas represent the summary statistic value
#' if (require(graphics)) {
#'   polygon(c(out1$x[1:400], 0.4, 0),
#'           c(out1$y[1:400], 1, 1),
#'           density = 10, angle = 60,
#'           border = "transparent")
#'   polygon(c(out1$x[401:999], 1, 0.4),
#'           c(out1$y[401:999], 0, 0),
#'           density = 10, angle = 60,
#'           border = "transparent")
#' }
#'
#' @export
conprof <- function(data, grpID = 1, n = 999, graph = TRUE, add = FALSE, ...) {
  if (inherits(data, "sf")) {
    data <- st_drop_geometry(data)
  }

  if (ncol(data) < 2 || !is.numeric(as.matrix(data))) {
    stop(
      "'data' must be a numeric matrix with at least two columns",
      call. = FALSE
    )
  }

  if (length(grpID) > 1) {
    warning(
      "'grpID' has more than one value, using the first value",
      call. = FALSE
    )
    grpID <- grpID[1]
  }

  colsum <- sum(data[, grpID])
  rowsum <- apply(data, 1, function(z) sum(z))

  xval <- rbind(seq(0, 1, length.out = n))
  yval <- numeric(n)
  threshold <- rowsum %*% xval

  for (i in 1:n) {
    INDEX <- (data[, grpID] >= threshold[, i])
    yval[i] <- sum(data[INDEX, grpID]) / colsum
  }

  val <- list("x" = c(xval, 1), "y" = c(yval, 0))

  if (graph) {
    if (!add) {
      plot(
        NA,
        xlim = c(0, 1.05),
        ylim = c(0, 1.05),
        xaxt = "n",
        yaxt = "n",
        ...
      )
      intrval <- seq(0, 1, by = 0.2)
      axistxt <- intrval * 100
      axis(side = 1, at = intrval, labels = axistxt, ...)
      axis(side = 2, at = intrval, labels = axistxt, ...)
    }
    lines(val$x, val$y, ...)
  }

  # Proportion of the group
  p <- sum(data[[grpID]]) / sum(data)

  above <- which(val$x >= p)
  below <- which(val$x < p)
  partA <- sum(val$y[above]) / n
  partB <- p - (sum(val$y[below]) / n)

  d <- (partB + partA) / (1 - p)

  return(list(x = val$x, y = val$y, d = d))
}

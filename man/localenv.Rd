\name{localenv}
\alias{localenv}
\title{Local Population Composition}
\description{
\code{localenv} calculates the local population composition at each data point from a matrix of coordinates, or an object of class \code{Spatial} or \code{ppp}.
}
\usage{
localenv(x, data, power = 3, 
         weighting = c("biweight", "unweighted", "inverse", "exponential"),
         normalize = TRUE, maxdist, sprel, useExp = NULL, scale = NULL, \dots)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{power}{a numeric value that determines the change rate of the biweight and inverse distance weighting functions. The default is 3. Typically 1-5.}
  \item{weighting}{a character string indicating the spatial weighting scheme. One of \dQuote{biweight} (default), \dQuote{unweighted}, \dQuote{inverse}, or \dQuote{exponential}.}
  \item{normalize}{logical. If TRUE, distances are divided by the bandwidth before weights are calculated so results are not affected by the coordinate measurement unit. The default is TRUE.}
  \item{maxdist}{an optional numeric value specifying a search radius for the construction of each local environment. Must be a positive value, or zero.}
  \item{sprel}{deprecated and no longer supported. Use coordinates or sf geometry instead.}
  \item{useExp}{deprecated and ignored. See \code{?localenv} for the current \sQuote{weighting} syntax.}
  \item{scale}{deprecated and ignored. See \code{?localenv} for the current \sQuote{normalize} syntax.}
  \item{\dots}{unused. The deprecated argument \code{tol} is accepted for compatibility and ignored.}
}
\details{
At each data point in \sQuote{x}, \code{localenv} calculates the weighted average of the populations of all points that are within a search radius \sQuote{maxdist}. The output from this function is an essential component to compute the spatial segregation measures.

If \sQuote{maxdist} is zero, no spatial weighting is applied: each local environment is set equal to the observed population values for the focal unit. This is useful when indices are desired without spatially weighted local environments.

The \sQuote{weighting} argument controls how neighbours within the bandwidth are weighted. The \dQuote{unweighted} option gives each neighbour equal weight. The \dQuote{biweight} option uses

\deqn{w(d) = (1 - (d / bw)^p)^p}{%
      w(d) = (1 - (d / bw)^p)^p}

where \emph{d} is the distance between the focal point and the neighbour, \emph{bw} is the bandwidth, and \emph{p} is \sQuote{power}. The \dQuote{inverse} option uses

\deqn{w(d^*) = \frac{1}{(d^* / bw)^p}}{%
      w(d*) = 1 / ((d* / bw)^p)}

where zero distances in \eqn{d} are replaced by half the minimum positive distance from focal point \eqn{i} to its neighbours within the bandwidth, and all nonzero distances are left unchanged. If there is no positive-distance neighbour within the bandwidth, the bandwidth itself is used for the zero-distance replacement. The \dQuote{exponential} option uses

\deqn{w(d) = e^{-d / bw}}{%
      w(d) = exp(-d / bw)}

By default, \emph{d} is normalized by the bandwidth before weights are calculated. If \sQuote{normalize} is FALSE, \emph{bw} is set to 1 in the weighting formulas. The bandwidth is still respected when selecting which neighbours are included.

If \sQuote{maxdist} is supplied, it is used as the bandwidth and as a search radius. If \sQuote{maxdist} is not supplied, the geometric mean of positive pairwise distances is used as the bandwidth. For large data sets, this bandwidth is estimated from a reproducible sample of points.

If \sQuote{maxdist} is not provided (default), the automatic bandwidth uses the geometric mean of pairwise distances, which is less sensitive to long-distance pairs than the arithmetic mean. To avoid the cost of a full pairwise distance matrix for large data sets, the automatic bandwidth is estimated from a reproducible sample of 500 points. It is still recommended to specify \sQuote{maxdist} when a substantively meaningful bandwidth is known, or to speed up the calculation process.

\sQuote{sprel} inputs of class \code{dist} or \code{nb} are no longer supported. Local environments are constructed from coordinates or sf geometry using the same C++ engine used by \code{\link{spseg}}.
}
\value{
An object of \code{\link{SegLocal-class}}.
}
\note{
Note that this function is not to interpolate between data points. The calculation of each local environment involves the point itself, so the simple inverse distance function with a power of 2 or more should be used with care.
}
\author{
Seong-Yun Hong
}
\seealso{
\code{\link{SegLocal-class}}, \code{\link{spseg}}, \code{dist}
}
\examples{
# uses the idealised landscapes in 'segdata'
data(segdata)
coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))

# complete segregation:
# default biweight kernel with normalized distance and power = 3
xx1 <- localenv(coords, data = segdata[,1:2])
print(xx1)

# inverse distance
xx2 <- localenv(coords, data = segdata[,1:2], weighting = "inverse")
print(xx2)

# unweighted local environments
xx3 <- localenv(coords, data = segdata[,1:2], weighting = "unweighted")
print(xx3)

# checkerboard pattern:
# biweight kernel with different p values
vv1 <- localenv(coords, data = segdata[,5:6], weighting = "biweight", power = 1)
vv2 <- localenv(coords, data = segdata[,5:6], weighting = "biweight", power = 2)
vv3 <- localenv(coords, data = segdata[,5:6], weighting = "biweight", power = 3)

# local environments can be exported as centroid coordinates plus values
head(as.data.frame(vv1))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_data <- rbind(
    transform(st_as_sf(vv1, columns = "C1"), power = "p = 1"),
    transform(st_as_sf(vv2, columns = "C1"), power = "p = 2"),
    transform(st_as_sf(vv3, columns = "C1"), power = "p = 3")
  )

  ggplot2::ggplot(plot_data) +
    ggplot2::geom_sf(ggplot2::aes(color = C1), size = 2) +
    ggplot2::facet_wrap(~ power) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_minimal()
}
}

\name{localenv}
\alias{localenv}
\alias{localenv.get}
\title{Local Population Composition}
\description{
\code{localenv} calculates the local population composition at each data point from a matrix of coordinates, or an object of class \code{Spatial} or \code{ppp}.
}
\usage{
localenv(x, data, power = 2, useExp = TRUE, scale = TRUE, maxdist, sprel, 
         tol = .Machine$double.eps)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{power}{a numeric value that determines the change rate of a distance weight function. If zero, all data points have the same weight regardless of the distance. Typically 1-5.}
  \item{useExp}{logical. If FALSE, use a simple inverse distance function instead of a negative exponential function. See \sQuote{Details}.}
  \item{scale}{logical. If TRUE, distances are divided by the bandwidth before weights are calculated so results are not affected by the coordinate measurement unit. The default is TRUE.}
  \item{maxdist}{an optional numeric value specifying a search radius for the construction of each local environment. Must be a positive value, or zero.}
  \item{sprel}{an optional object of class \code{dist} or \code{nb}. See \sQuote{Details}.}
  \item{tol}{a small, positive non-zero value. If \sQuote{useExp} is FALSE, this value will be added to the denominator to avoid the divide-by-zero error.}
}
\details{
At each data point in \sQuote{x}, \code{localenv} calculates the weighted average of the populations of all points that are within a search radius \sQuote{maxdist}. The output from this function is an essential component to compute the spatial segregation measures.

By default, the weight of each point is calculated from a negative exponential function, which is defined as:

\deqn{w(d) = e^{-d \times power}}{%
      w(d) = exp(-d * power)}

where \emph{d} is the Euclidean distance between two points. By default, \emph{d} is normalized by the bandwidth before weights are calculated. If \sQuote{maxdist} is supplied, it is used as the bandwidth and as a search radius. If \sQuote{maxdist} is not supplied, the geometric mean of positive pairwise distances is used as the bandwidth. For large data sets, this bandwidth is estimated from a reproducible sample of points.

If \sQuote{useExp} is FALSE, a simple inverse distance function is used to calculate the weight of each point:

\deqn{w(d) = \frac{1}{(d + error)^{power}}}{%
      w(d) = 1 / (d + error)^power}

If \sQuote{maxdist} is not provided (default), the automatic bandwidth uses the geometric mean of pairwise distances, which is less sensitive to long-distance pairs than the arithmetic mean. To avoid the cost of a full pairwise distance matrix for large data sets, the automatic bandwidth is estimated from a reproducible sample of 500 points. It is still recommended to specify \sQuote{maxdist} when a substantively meaningful bandwidth is known, or to speed up the calculation process.

If a distance measure other than the Euclidean distance is required to represent spatial proximity between the points, users can provide an object of class \code{dist}, which contains the distances between all pairs of the points, through an optional argument \sQuote{sprel}. One convenient way of obtaining such information may be the use of the function \code{dist}, which offers a variety of distance measures, such as Manhattan, Canberra, and Minkowski.

Or alternatively, one can supply an object of class \code{nb} to use a k-nearest neighbour averaging or polygon contiguity.
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
\code{\link{SegLocal-class}}, \code{\link{spatseg}}, \code{dist}
}
\examples{
# uses the idealised landscapes in 'segdata'
data(segdata)
coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))

# complete segregation:
# negative exponential function of distance
xx1 <- localenv(coords, data = segdata[,1:2])
print(xx1)

# inverse distance
xx2 <- localenv(coords, data = segdata[,1:2], useExp = FALSE)
print(xx2)

# inverse distance with p = 0, i.e., weight all points equally
xx3 <- localenv(coords, data = segdata[,1:2], useExp = FALSE, power = 0)
print(xx3)

\dontrun{
# checkerboard pattern:
# negative exponential function with different p values
vv1 <- localenv(coords, data = segdata[,5:6], power = 1)
spplot(vv1, main = "Negative exponetial with p = 1")

vv2 <- localenv(coords, data = segdata[,5:6])
spplot(vv2, main = "Negative exponetial with p = 2")

vv3 <- localenv(coords, data = segdata[,5:6], power = 3)
spplot(vv3, main = "Negative exponetial with p = 3")}
}

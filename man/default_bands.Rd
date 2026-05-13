\name{default_bands}
\alias{default_bands}
\title{Default Distance Bands}
\description{
Calculates default distance bands from sampled pairwise distances.
}
\usage{
default_bands(x, data = NULL, n = 500, probs = seq(0.1, 0.9, 0.1),
              weighted = TRUE)
}
\arguments{
  \item{x}{a numeric coordinate matrix or a spatial object accepted by \code{\link{localenv}}.}
  \item{data}{optional population data. If omitted and \code{x} contains attached data, those data are used.}
  \item{n}{maximum number of points to sample.}
  \item{probs}{probabilities passed to \code{\link{quantile}}.}
  \item{weighted}{logical. If TRUE, sample points using row sums of \code{data} as population weights.}
}
\details{
This helper samples up to \code{n} points, calculates pairwise distances among the sampled points, and returns the requested quantiles. It is used by \code{\link{spseg}} when non-legacy output is requested and no bandwidths are supplied.
}
\value{
A numeric vector of distance bands.
}
\seealso{
\code{\link{spseg}}, \code{\link{localenv}}
}

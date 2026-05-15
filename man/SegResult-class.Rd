\name{SegResult-class}
\docType{class}
\alias{SegResult-class}
\alias{SegResult}
\alias{show,SegResult-method}
\alias{print.SegResult}
\alias{as.list.SegResult}
\alias{as.data.frame.SegResult}
\title{Unified Segregation Result Class}
\description{
Container for local environments, spatial segregation indices, or both, over one or more bandwidths.
}
\usage{
SegResult(coords = NULL, data = NULL, env = NULL, bands, indices = list(),
          measures = character(), weighting = character(), power = numeric(),
          normalize = logical(), neighbors = list(type = "radius"),
          proj4string = st_crs(as.character(NA)), call = NULL)
}
\arguments{
  \item{coords}{a numeric coordinate matrix, or NULL when coordinates are not stored.}
  \item{data}{a numeric population matrix, or NULL when input data are not stored.}
  \item{env}{a list of local-environment matrices, one per band, or NULL.}
  \item{bands}{numeric vector of bandwidths/search radii.}
  \item{indices}{a list of segregation index results. Overall measures are stored in \code{indices$overall}; pairwise measures, when requested, are stored in \code{indices$pairwise}.}
  \item{measures}{character vector of calculated measures.}
  \item{weighting}{character string naming the spatial weighting scheme.}
  \item{power}{numeric weighting power.}
  \item{normalize}{logical; whether distances were normalized by bandwidth.}
  \item{neighbors}{a list describing the neighborhood strategy. Currently radius bands are supported.}
  \item{proj4string}{coordinate reference system.}
  \item{call}{the matched call that produced the object, or NULL.}
}
\details{
\code{SegResult} is the unified result object used by \code{\link{spseg}} for multi-band calculations, pairwise comparisons, and non-legacy output modes. It can represent local environments only, indices only, or both. When \code{output = "indices"} in \code{\link{spseg}}, coordinates, input data, and local environments are not stored.

\code{as.data.frame(x, what = "indices")} returns long-form overall scalar measures. \code{as.data.frame(x, what = "pairwise")} returns long-form pairwise scalar measures.
}
\value{
An object of class \code{SegResult}.
}
\seealso{
\code{\link{spseg}}, \code{\link{SegSpatial-class}}, \code{\link{SegLocal-class}}
}

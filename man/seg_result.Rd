\name{seg_result}
\alias{seg_result}
\alias{print.seg_result}
\alias{as.list.seg_result}
\alias{as.data.frame.seg_result}
\alias{st_as_sf.seg_result}
\title{Spatial Segregation Result Objects}
\description{
S3 result objects returned by \code{\link{spseg}}.
}
\usage{
\method{print}{seg_result}(x, \dots)
\method{as.list}{seg_result}(x, \dots)
\method{as.data.frame}{seg_result}(x, row.names = NULL, optional = FALSE,
  \dots, what = c("indices", "pairwise", "exposure", "env"),
  bands = NULL, columns = NULL, shape = c("wide", "long"))
\method{st_as_sf}{seg_result}(x, \dots, bands = NULL, columns = NULL,
  shape = c("wide", "long"))
}
\arguments{
  \item{x}{A \code{seg_result} object.}
  \item{row.names, optional}{Arguments passed from the
  \code{as.data.frame()} generic.}
  \item{\dots}{Ignored.}
  \item{what}{Which part of the result to extract.}
  \item{bands}{Optional subset of bands to extract. Values must match
  \code{x$bands}.}
  \item{columns}{Optional subset of local-environment columns to extract.}
  \item{shape}{Whether local environments should be returned in wide or long
  format.}
}
\details{
A \code{seg_result} is a plain S3 list. The most commonly used fields are:
\code{bands}, \code{indices}, \code{env}, \code{coords}, \code{data},
\code{geometry}, \code{neighbors}, and \code{crs}. When
\code{spseg(output = "indices")}, local environments and input data are not
stored. When \code{spseg(output = "localenv")}, indices are not calculated.
}
\value{
\code{as.data.frame()} returns a data frame. \code{sf::st_as_sf()} returns an
\code{sf} object containing stored local environments.
}
\seealso{
\code{\link{spseg}}
}

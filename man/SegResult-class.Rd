\name{SegResult-class}
\docType{class}
\alias{SegResult-class}
\alias{SegResult}
\alias{show,SegResult-method}
\alias{print.SegResult}
\alias{as.list.SegResult}
\alias{as.data.frame.SegResult}
\alias{st_as_sf.SegResult}
\title{Unified Segregation Result Class}
\description{
Container for local environments, spatial segregation indices, or both, over one or more bandwidths.
}
\usage{
SegResult(coords = NULL, data = NULL, env = NULL, bands, indices = list(),
          measures = character(), weighting = character(), power = numeric(),
          normalize = logical(), neighbors = list(type = "radius"),
          geometry = NULL, proj4string = st_crs(as.character(NA)),
          call = NULL)
}
\arguments{
  \item{coords}{a numeric coordinate matrix, or NULL when coordinates are not stored.}
  \item{data}{a numeric population matrix, or NULL when input data are not stored.}
  \item{env}{a list of local-environment matrices, one per band, or NULL.}
  \item{geometry}{an optional \code{sfc} geometry object with one feature per row.}
  \item{bands}{numeric vector of bandwidths/search radii.}
  \item{indices}{a list of segregation index results. Overall measures are stored in \code{indices$overall}; pairwise measures, when requested, are stored in \code{indices$pairwise}.}
  \item{measures}{character vector of calculated measures.}
  \item{weighting}{character string naming the spatial weighting scheme.}
  \item{power}{numeric weighting power.}
  \item{normalize}{logical; whether distances were normalized by bandwidth.}
  \item{neighbors}{a list describing the neighborhood strategy, threshold values, threshold units, and search engine. Radius neighborhoods use distance thresholds; kNN neighborhoods use population-count thresholds.}
  \item{proj4string}{coordinate reference system.}
  \item{call}{the matched call that produced the object, or NULL.}
}
\details{
\code{SegResult} is the unified result object used by \code{\link{spseg}} for multi-band calculations, pairwise comparisons, and non-legacy output modes. It can represent local environments only, indices only, or both. When \code{output = "indices"} in \code{\link{spseg}}, coordinates, input data, and local environments are not stored.

\code{as.data.frame(x, what = "indices")} returns long-form overall scalar measures. \code{as.data.frame(x, what = "pairwise")} returns long-form pairwise scalar measures. \code{as.data.frame(x, what = "env")} returns centroid coordinates and stored local environment values. Use \code{bands}, \code{columns}, and \code{format = "wide"} or \code{format = "long"} to select and shape the environment output.

\code{sf::st_as_sf(x)} converts stored local environments to an \code{sf} object. Use \code{bands} to select bandwidths, \code{columns} to select local environment columns, and \code{format = "wide"} or \code{format = "long"} to choose one row per feature or one row per feature-band combination.
}
\value{
An object of class \code{SegResult}.
}
\seealso{
\code{\link{spseg}}, \code{\link{SegSpatial-class}}, \code{\link{SegLocal-class}}
}
\examples{
data(segdata)
coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))

result <- spseg(coords, data = segdata[, 1:2], bands = c(1, 2),
                output = "full")

# Overall scalar indices by bandwidth
as.data.frame(result, what = "indices")

# Stored local environments as a coordinate data frame
head(as.data.frame(result, what = "env", columns = "A1"))

# Long form is convenient for plotting or grouped summaries
head(as.data.frame(result, what = "env", columns = "A1",
                   format = "long"))

# Stored local environments as sf point geometry
if (requireNamespace("ggplot2", quietly = TRUE)) {
  result_sf <- st_as_sf(result, bands = 1, columns = "A1")
  ggplot2::ggplot(result_sf) +
    ggplot2::geom_sf(ggplot2::aes(color = bw_1_A1), size = 2) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_minimal()
}
}

\name{spseg}
\alias{spseg}
\title{Spatial Segregation Measures}
\description{
Calculates spatial segregation measures from point coordinates or \code{sf}
geometries using local population environments.
}
\usage{
spseg(x, data = NULL, measures = "all", bands = NULL,
      weighting = c("biweight", "unweighted", "inverse", "exponential"),
      power = 3, normalize = TRUE, neighbors = c("radius", "knn"),
      search = c("kdtree", "brute"),
      surface = c("raw", "grid", "pycno"),
      comparison = c("overall", "pairwise", "both"),
      output = c("indices", "full", "localenv"),
      verbose = FALSE, \dots)
}
\arguments{
  \item{x}{A numeric matrix or data frame with two coordinate columns, or an
  object of class \code{sf}. For polygon \code{sf} inputs, centroids are used as
  the focal coordinates.}
  \item{data}{A numeric matrix or data frame with one row per coordinate or
  feature and at least two population-group columns. If \code{x} is an
  \code{sf} object, \code{data} may be omitted and the non-geometry columns of
  \code{x} are used.}
  \item{measures}{A character vector naming measures to compute: \dQuote{all},
  \dQuote{exposure}, \dQuote{information}, \dQuote{diversity}, and
  \dQuote{dissimilarity}.}
  \item{bands}{Neighborhood thresholds. With \code{neighbors = "radius"}, these
  are distance bandwidths in the units of \code{x}. With
  \code{neighbors = "knn"}, these are population-count thresholds. If omitted
  for radius neighborhoods, \code{\link{default_bands}} is used.}
  \item{weighting}{Spatial weighting scheme used within each neighborhood. One
  of \dQuote{biweight}, \dQuote{unweighted}, \dQuote{inverse}, or
  \dQuote{exponential}.}
  \item{power}{Power parameter for the biweight and inverse-distance weighting
  schemes.}
  \item{normalize}{Logical. If \code{TRUE}, distances are divided by the
  bandwidth before weights are calculated. The bandwidth still defines the
  search radius when \code{neighbors = "radius"}.}
  \item{neighbors}{Neighborhood strategy. \dQuote{radius} uses distance
  bandwidths. \dQuote{knn} uses count-based adaptive neighborhoods.}
  \item{search}{Neighbor-search backend. \dQuote{kdtree} uses the exact
  nanoflann kd-tree backend. \dQuote{brute} computes candidate distances
  directly and is mainly useful for validation.}
  \item{surface}{Population surface used before local environments are
  calculated. \dQuote{raw} uses the input support directly. \dQuote{grid}
  redistributes areal counts over a square grid. \dQuote{pycno} applies
  groupwise pycnophylactic smoothing on that grid while preserving source
  polygon group totals.}
  \item{comparison}{Whether to calculate \dQuote{overall} multi-group measures,
  \dQuote{pairwise} two-group measures, or \dQuote{both}.}
  \item{output}{Return mode. \dQuote{indices} stores only segregation indices.
  \dQuote{full} stores both indices and local environments. \dQuote{localenv}
  stores local environments without calculating indices.}
  \item{verbose}{Logical. If \code{TRUE}, report progress messages.}
  \item{\dots}{Additional population-surface options, such as \code{cellsize},
  \code{celldim}, \code{nrow}, \code{ncol}, \code{max_iter}, or
  \code{converge}.}
}
\details{
\code{spseg()} is the primary interface for the package. Local environments
and segregation indices are produced by a unified C++ engine. The default
calculation uses normalized biweight local environments with \code{power = 3}
and radius neighborhoods.

By default, \code{spseg()} returns indices over default distance bands. These
default bands are the first through fifth deciles of distances among a
reproducible population-weighted sample of points; see
\code{\link{default_bands}}. If \code{output = "indices"}, coordinates, input
data, and local environments are not stored in the result.

Set \code{bands = 0} to compute indices without spatial weighting. In that
case, each local environment is a copy of the observed population values for
the focal unit.

When \code{neighbors = "knn"}, \code{bands} is interpreted as population-count
thresholds. The focal unit is included fully. If the focal unit's population
exceeds a threshold, its observed composition is used for that threshold.
Otherwise, neighboring units are added in order of centroid distance until the
threshold is reached, with the final boundary unit included fractionally when
needed.

Overall measures are stored in \code{result$indices$overall}. Pairwise
measures, when requested, are stored in \code{result$indices$pairwise}.
Pairwise scalar measures are returned as one symmetric matrix per bandwidth
with one row and column for each population group. Units with no population in
either group for a given pair do not contribute to that pair's accumulators.

The \code{surface} argument controls population-surface construction before
local environments are calculated. \code{surface = "grid"} redistributes
areal counts over a square grid. \code{surface = "pycno"} applies internal
groupwise pycnophylactic interpolation: each group is smoothed separately by
averaging each grid cell with neighboring cells, then rescaled within each
source polygon after every iteration so known group totals are preserved.

When computing the spatial information theory index, zero proportions are
evaluated using the limiting value \eqn{x \log(x) = 0} as \eqn{x} approaches
zero. This avoids artificial perturbation of the data and prevents
\code{0 * log(0)} from producing \code{NaN}.
}
\value{
An S3 object of class \code{seg_result}. It is a list with fields including
\code{bands}, \code{indices}, \code{env}, \code{coords}, \code{data},
\code{geometry}, \code{neighbors}, and \code{crs}. Some fields are
\code{NULL} depending on \code{output}.
}
\note{
The exposure/isolation index, P, is presented in matrix form. The spatial
exposure of group \sQuote{m} to group \sQuote{n} is located in row
\sQuote{m} and column \sQuote{n}. The matrix is rarely symmetric in practice.
}
\author{
Seong-Yun Hong
}
\references{
Reardon, S. F. and O'Sullivan, D. (2004) Measures of spatial segregation.
\emph{Sociological Methodology}, \bold{34}, 121-162.

Reardon, S. F., Farrell, C. R., Matthews, S. A., O'Sullivan, D., Bischoff, K.,
and Firebaugh, G. (2009) Race and space in the 1990s: Changes in the
geographic scale of racial residential segregation, 1990-2000.
\emph{Social Science Research}, \bold{38}, 55-70.
}
\seealso{
\code{\link{default_bands}}, \code{\link{as.data.frame.seg_result}}
}
\examples{
data(segdata)

coords <- as.matrix(expand.grid(x = seq(0.5, 9.5), y = seq(0.5, 9.5)))
test_data <- segdata[, 1:2]

# Default behavior: indices across default bands.
idx <- spseg(coords, data = test_data)
idx
as.data.frame(idx)

# Store local environments for selected distance bands.
env <- spseg(coords, data = test_data, bands = c(0, 2), output = "localenv")
head(as.data.frame(env, what = "env", shape = "long"))

# Store indices and local environments together.
full <- spseg(coords, data = test_data, bands = c(1, 2), output = "full")
as.data.frame(full, what = "exposure")

# sf inputs preserve geometry in local-environment outputs.
geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)))) |>
  sf::st_make_grid(cellsize = 1)
grd_sf <- sf::st_sf(test_data, geometry = geom)

sf_env <- spseg(grd_sf, bands = 2, output = "localenv")
sf::st_as_sf(sf_env)
}

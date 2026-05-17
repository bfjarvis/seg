\name{spseg}
\alias{spseg}
\alias{spatseg}
\title{Spatial Segregation Measures}
\description{
Calculates the set of spatial segregation measures developed by Reardon and O'Sullivan (2004).
}
\usage{
spseg(x, data, measures = "all", sprel = NULL, maxdist = NULL,
      bands = NULL,
      weighting = c("biweight", "unweighted", "inverse", "exponential"),
      power = 3, normalize = TRUE, neighbors = c("radius", "knn"),
      search = c("kdtree", "brute"),
      smoothing = NULL,
      comparison = c("overall", "pairwise", "both"),
      output = c("indices", "full", "localenv", "legacy"),
      useC = TRUE, negative.rm = FALSE, 
      verbose = FALSE, \dots)
spatseg(env, method = "all", useC = TRUE, negative.rm = FALSE, \dots)
}
\arguments{
  \item{x}{an object of \code{\link{SegLocal-class}}, a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{env}{an object of \code{\link{SegLocal-class}}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (i.e., at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{measures}{a vector of character strings indicating the measure or measures to be computed. This must be one or more of the strings \dQuote{all} (default), \dQuote{exposure}, \dQuote{information}, \dQuote{diversity}, and \dQuote{dissimilarity}. Abbreviations are accepted, as long as it is clear which measure is meant.}
  \item{method}{deprecated synonym for \code{measures}, retained for compatibility in \code{spatseg} and accepted through \dots in \code{spseg}.}
  \item{sprel}{deprecated and no longer supported. Use coordinates or sf geometry instead.}
  \item{maxdist}{an optional numeric value passed to \code{\link{localenv}} specifying the bandwidth/search radius for local environments. If zero, no spatial weighting is applied and the observed data are used as local environments.}
  \item{bands}{an optional numeric vector of neighborhood thresholds. With \code{neighbors = "radius"}, these are distance bandwidths/search radii and are a clearer synonym for supplying multiple values through \code{maxdist}. If both \code{bands} and \code{maxdist} are NULL, default radius bands are computed by \code{\link{default_bands}}, except when \code{output = "legacy"} is explicitly requested. With \code{neighbors = "knn"}, \code{bands} must be supplied and gives population-count thresholds.}
  \item{weighting}{a character string passed to \code{\link{localenv}} indicating the spatial weighting scheme. One of \dQuote{biweight} (default), \dQuote{unweighted}, \dQuote{inverse}, or \dQuote{exponential}.}
  \item{power}{a numeric value passed to \code{\link{localenv}} that determines the change rate of the biweight and inverse distance weighting functions. The default is 3.}
  \item{normalize}{logical. If TRUE, distances are divided by the bandwidth before local environment weights are calculated. The default is TRUE.}
  \item{neighbors}{a character string indicating the neighborhood strategy. \dQuote{radius} uses distance bandwidths/search radii and is the default. \dQuote{knn} uses count-based adaptive neighborhoods; in this mode \code{bands} gives population-count thresholds rather than distances.}
  \item{search}{a character string indicating the neighbor-search backend used by the unified engine. \dQuote{kdtree} uses the nanoflann kd-tree backend and is the default. \dQuote{brute} computes all centroid distances directly, which is useful for validation and can be competitive when neighborhoods are very dense.}
  \item{smoothing}{an optional named list controlling population smoothing before local environments are constructed. If NULL, smoothing is not used. See \sQuote{Details}.}
  \item{comparison}{a character string indicating whether to calculate \dQuote{overall} measures using all supplied groups, \dQuote{pairwise} measures for every two-group comparison, or \dQuote{both}.}
  \item{output}{a character string controlling the returned object. \dQuote{indices} is the default and returns only segregation indices across the chosen or default bands without storing coordinates, input data, or local environments. \dQuote{full} returns a \code{\link{SegResult-class}} with both local environments and indices. \dQuote{localenv} stores local environments but does not calculate indices. \dQuote{legacy} preserves the previous scalar-band \code{\link{SegSpatial-class}} return when explicitly requested.}
  \item{useC}{deprecated logical retained for compatibility. Segregation calculations now use the C++ engine regardless of this value.}
  \item{negative.rm}{logical. If TRUE, all geographic units where at least one group (i.e., column) has a negative population will be removed. If FALSE, negative values will be replaced with zero.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
  \item{\dots}{additional arguments passed to \code{\link{localenv}} when \sQuote{x} is raw spatial data. Deprecated top-level smoothing arguments \code{nrow}, \code{ncol}, \code{nocol}, \code{window}, and \code{sigma} are accepted for compatibility with a warning and folded into \sQuote{smoothing}. The deprecated arguments \code{method} and \code{tol} are accepted for compatibility; \code{method} is used as \code{measures} when \code{measures} is not supplied, and \code{tol} is ignored.}
}
\details{
\code{spseg} computes the set of spatial segregation measures proposed by Reardon and O'Sullivan. It can be called directly on a \code{\link{SegLocal-class}} object, or on raw spatial data from which \code{spseg} first constructs local environments.

\code{spatseg} is deprecated and retained as a compatibility wrapper for existing code. New code should call \code{spseg} instead.

Local environment construction is the main spatial modelling step. When \code{spseg} is called on raw spatial data, it first constructs local environments. By default, \code{spseg} uses normalized biweight local environments with \code{power = 3}, matching the defaults in \code{\link{localenv}}. The most important local environment choices are \code{maxdist} or \code{bands} (the bandwidth/search radius), \code{weighting}, \code{power}, and \code{normalize}. See \code{\link{localenv}} for details on how these choices define local population composition.

By default, \code{spseg} returns a \code{\link{SegResult-class}} object containing indices across default distance bands. The default bands are the first through fifth deciles of distances among a reproducible population-weighted sample of points; see \code{\link{default_bands}}. Local environments and segregation indices are produced by a unified C++ engine. If \code{output = "indices"}, local environments are used only transiently inside the calculation and are not stored in the result.

Set \code{output = "legacy"} to request the previous scalar-band \code{\link{SegSpatial-class}} return. If neither \code{maxdist} nor \code{bands} is supplied in legacy mode, a single localenv-style geometric-mean bandwidth is used.

Overall measures are stored in \code{indices$overall}. Pairwise measures, when requested, are stored in \code{indices$pairwise}. Pairwise scalar measures are returned as one symmetric matrix per bandwidth with one row and column for each population group. For pairwise calculations, each two-group comparison is evaluated within its own two-group population universe; units with no population in either group for a given pair do not contribute to that pair's accumulators. If either group has zero total population over the whole study area, the corresponding pairwise value is \code{NA}.

Set \code{maxdist = 0} to compute indices without spatially weighted local environments. In this case, each local environment is a copy of the observed data for the focal unit.

When \code{neighbors = "knn"}, \code{bands} is interpreted as population-count thresholds. For each focal unit, the focal unit is included fully. If the focal unit's population exceeds a threshold, its observed composition is used for that threshold. Otherwise, neighbouring units are added in order of centroid distance until the threshold is reached, with the final boundary unit included fractionally when necessary. This treats each areal unit as representing the local environment of the people in that unit, without assuming an internal population distribution within the geometry.

The \code{search} argument controls how candidate neighbours are found in the unified engine. The default \code{search = "kdtree"} uses a kd-tree for exact neighbor lookup. \code{search = "brute"} uses direct all-pairs distance evaluation and is mainly intended as a transparent validation backend or for cases where kd-tree overhead is not worthwhile.

Population smoothing before local environment construction is also possible through the \sQuote{smoothing} argument, but this is a secondary option. If supplied, \sQuote{smoothing} should be a named list with entries corresponding to the older smoothing arguments, such as \code{smoothing = list(smoothing = "equal", nrow = 100, ncol = 100)} or \code{smoothing = list(smoothing = "kernel", nrow = 100, ncol = 100, window = window, sigma = sigma)}. Consult \code{\link{localenv}} first when choosing the spatial scale and weighting scheme for the local environment calculation.

When computing the spatial information theory index, zero proportions are evaluated using the limiting value \eqn{x \log(x) = 0} as \eqn{x} approaches zero. This avoids artificial perturbation of the data and prevents \code{0 * log(0)} from producing \code{NaN}.

The spatial dissimilarity index is calculated using the generalized local-environment formula from Reardon and O'Sullivan. For two groups, this is algebraically equivalent to the pairwise Duncan-Duncan expression applied to local environment proportions. With more than two groups, it is a multi-group generalization rather than a first-two-groups pairwise measure.
}
\value{
An object of \code{\link{SegResult-class}} by default, or \code{\link{SegSpatial-class}} when \code{output = "legacy"} is explicitly requested.
}
\note{
For entropy calculations, terms of the form \eqn{x \log(x)} are evaluated as zero when \eqn{x = 0}, following the limiting value as \eqn{x} approaches zero. The deprecated \code{tol} argument is therefore ignored when supplied.

The exposure/isolation index, P, is presented in a matrix form. The spatial exposure of group \sQuote{m} to group \sQuote{n} is located in the row \sQuote{m} and column \sQuote{n} of the matrix. The matrix is rarely symmetric in practice so the spatial exposure index should be interpreted with care.

The spatial isolation index values are given in the diagonal cells of the matrix; cell value at (m, m) indicates the degree of spatial isolation for group \sQuote{m} for example.
}
\author{
Seong-Yun Hong
}
\references{
Reardon, S. F. and O'Sullivan, D. (2004) Measures of spatial segregation. \emph{Sociological Methodology}, \bold{34}, 121-162.

Reardon, S. F., Farrell, C. R., Matthews, S. A., O'Sullivan, D., Bischoff, K., and Firebaugh, G. (2009) Race and space in the 1990s: Changes in the geographic scale of racial residential segregation, 1990-2000. \emph{Social Science Research}, \bold{38}, 55-70.
}
\seealso{
\code{\link{SegSpatial-class}}, \code{\link{localenv}}, \code{kernel2d}
}
\examples{
# uses the idealised landscapes in 'segdata'
data(segdata)
geom <- st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)))) |>
  st_make_grid(cellsize = 1)
grd_sf <- st_sf(geometry = geom)
test.df <- segdata[,1:2]

# explicit workflow: compute local environments, then spatial indices
env <- localenv(grd_sf, data = test.df, maxdist = 2)
xx1 <- spseg(env)
print(xx1, digits = 3)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  env_sf <- st_as_sf(env, columns = "A1")
  ggplot2::ggplot(env_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = A1), color = "white") +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_void()
}

# selected measures can be requested by name
xx2 <- spseg(env, measures = c("information", "dissimilarity"))
as.list(xx2)

# one-call workflow using a data frame of coordinates
coords_df <- data.frame(x = rep(seq(0.5, 9.5), each = 10),
                        y = rep(seq(0.5, 9.5), times = 10))
xx3 <- spseg(
  coords_df,
  data = test.df,
  measures = c("information", "dissimilarity"),
  maxdist = 2,
  weighting = "biweight",
  power = 3,
  normalize = TRUE
)
as.list(xx3)

# optional smoothing is supplied as a named list
xx4 <- spseg(
  grd_sf,
  data = test.df,
  smoothing = list(smoothing = "equal", nrow = 10, ncol = 10),
  maxdist = 2
)
as.list(xx4)

# multiple bandwidths return a SegResult
xx5 <- spseg(
  coords_df,
  data = test.df,
  bands = c(1, 2, 3),
  output = "full"
)
as.data.frame(xx5)
head(as.data.frame(xx5, what = "env", columns = "A1"))
head(as.data.frame(xx5, what = "env", columns = "A1", shape = "long"))

xx6 <- spseg(
  coords_df,
  data = segdata[, 1:3],
  bands = c(1, 2),
  comparison = "both",
  output = "indices"
)
as.data.frame(xx6, what = "pairwise")

# count-based k-nearest-neighbour local environments
xx7 <- spseg(
  coords_df,
  data = test.df,
  bands = c(10, 25, 50),
  neighbors = "knn",
  search = "kdtree",
  weighting = "unweighted",
  output = "full"
)
head(as.data.frame(xx7, what = "env", shape = "long"))
}

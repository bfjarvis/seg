\name{spseg}
\alias{spseg}
\alias{spatseg}
\title{Spatial Segregation Measures}
\description{
Calculates the set of spatial segregation measures developed by Reardon and O'Sullivan (2004).
}
\usage{
spseg(x, data, method = "all", smoothing = "none", nrow = 100, ncol = 100, 
      window, sigma, useC = TRUE, negative.rm = FALSE, 
      verbose = FALSE, \dots)
spatseg(env, method = "all", useC = TRUE, negative.rm = FALSE, \dots)
}
\arguments{
  \item{x}{an object of \code{\link{SegLocal-class}}, a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{env}{an object of \code{\link{SegLocal-class}}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (i.e., at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{method}{a vector of character strings indicating an measure or measures to be computed. This must be one or more of the strings \dQuote{all} (default), \dQuote{exposure}, \dQuote{information}, \dQuote{diversity}, and \dQuote{dissimilarity}. Abbreviations are accepted, as long as it is clear which method is meant.}
  \item{smoothing}{a character string indicating how to perform spatial smoothing of the population data. This must be (an abbreviation of) one of the strings \dQuote{none} (default), \dQuote{kernel}, or \dQuote{equal}.}
  \item{nrow}{a numeric value indicating the number of row cells in the rasterised data surface. Ignored if \sQuote{smoothing} is \dQuote{none}.}
  \item{ncol}{a numeric value indicating the number of column cells.}
  \item{window}{an optional object of class \code{matrix} to be passed to \code{kernel2d}. See \sQuote{Details} in \code{\link{deseg}}.}
  \item{sigma}{an optional numeric value specifying the kernel bandwidth to be passed to \code{kernel2d}. See also \sQuote{Details} in \code{\link{deseg}}.}
  \item{useC}{logical. If TRUE, calculate the segregation values in C.}
  \item{negative.rm}{logical. If TRUE, all geographic units where at least one group (i.e., column) has a negative population will be removed. If FALSE, negative values will be replaced with zero.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
  \item{\dots}{optional arguments to be passed to \code{\link{localenv}} to compute the population composition of each local environment. The deprecated argument \code{tol} is accepted for compatibility and ignored.}
}
\details{
\code{spseg} computes the set of spatial segregation measures proposed by Reardon and O'Sullivan. It can be called directly on a \code{\link{SegLocal-class}} object, or on spatial data from which \code{spseg} first constructs local environments.

\code{spatseg} is deprecated and retained as a compatibility wrapper for existing code. New code should call \code{spseg} instead.

When \code{spseg} is called on raw spatial data, the population density surface is constructed by assuming that the population density is uniform in each census tract, or by using \code{kernel2d} in the package \pkg{splancs}. In the previous version (< 0.4.1), the function \code{rasterize} in \pkg{raster} was used in the case of the former, and \code{density.ppp} in \pkg{spatstat} for the latter.

When computing the spatial information theory index, zero proportions are evaluated using the limiting value \eqn{x \log(x) = 0} as \eqn{x} approaches zero. This avoids artificial perturbation of the data and prevents \code{0 * log(0)} from producing \code{NaN}.
}
\value{
An object of \code{\link{SegSpatial-class}}.
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

# no spatial smoothing: compute local environments and spatial indices
env <- localenv(grd_sf, data = test.df, maxdist = 2)
xx1 <- spseg(env)
print(xx1, digits = 3)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  env_sf <- cbind(grd_sf, as.data.frame(env@env))
  ggplot2::ggplot(env_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = A1), color = "white") +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_void()
}

# selected measures can be requested by name
xx2 <- spseg(env, method = c("information", "dissimilarity"))
as.list(xx2)
}

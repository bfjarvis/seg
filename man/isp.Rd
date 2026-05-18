\name{isp}
\alias{isp}
\title{Index of Spatial Proximity}
\description{
Computes the index of spatial proximity developed by White (1983). This measure estimates the level of clustering by comparing the average distance between members of the same group with that between all individuals (regardless of the groups to which they belong). The results may change drastically depending on the definition of distance.
}
\usage{
isp(x, data, nb, fun, verbose = FALSE, \dots)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{sf}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of geographic units in \sQuote{x}, and the number of columns should be greater than one (i.e., at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{nb}{an optional \code{matrix} object indicating the distances between the geographic units.}
  \item{fun}{a function for the calculation of proximity. The function should take a numeric vector as an argument (distance) and return a vector of the same length (proximity). If this is not specified, a negative exponential function is used by default.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
  \item{\dots}{optional arguments to be passed to \code{dist} when calculating the distances between the geographic units in \sQuote{x}. Ignored if \sQuote{nb} is given. See \code{help(dist)} for available options.}
}
\details{
\sQuote{nb} must be a square matrix (i.e., one that has the same number of rows and columns) but does not have to be symmetric. If \sQuote{nb} is not given, \code{whiteseg} attempts to create a distance matrix of \sQuote{x} using the function \code{dist} in \pkg{stats} and use it as \sQuote{nb}. The optional arguments in \sQuote{\dots} will be passed to \code{dist}.
}
\value{
A single numeric value indicating the degree of segregation; a value of 1 indicates absence of segregation, and values greater than 1.0 indicate clustering. If the index value is less than one, it indicates an unusual form of segregation (i.e., people live closer to other population groups).
}
\author{
Seong-Yun Hong
}
\references{
White, M. J. (1983). The measurement of spatial segregation. \emph{The American Journal of Sociology}, \bold{88}, 1008-1018.
}
\seealso{
\code{\link{dissim}}, \code{dist}
}
\examples{
# uses the idealised landscapes in 'segdata'
data(segdata)
geom <- sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,0)))) |> 
  sf::st_make_grid(cellsize = 1)
grd_sf <- sf::st_sf(geom)

d <- rep(NA, 8) # index of dissimilarity
p <- rep(NA, 8) # index of spatial proximity
for (i in 1:8) {
  idx <- 2 * i
  d[i] <- dissim(data = segdata[,(idx-1):idx])$d
  p[i] <- isp(grd_sf, data = segdata[,(idx-1):idx])
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_data <- do.call(
    rbind,
    lapply(1:8, function(i) {
      idx <- 2 * i - 1
      data.frame(
        pattern = paste0(
          LETTERS[i], ": D = ", round(d[i], 2), ", P = ", round(p[i], 2)
        ),
        group_share = segdata[, idx],
        geometry = sf::st_geometry(grd_sf)
      )
    })
  )
  plot_data <- sf::st_sf(plot_data)

  ggplot2::ggplot(plot_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = group_share), color = "white") +
    ggplot2::facet_wrap(~ pattern) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_void()
}
}

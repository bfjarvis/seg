\name{SegLocal-class}
\docType{class}
\alias{SegLocal-class}
\alias{coerce,SegLocal,SpatialPoints-method}
\alias{coerce,SegLocal,SpatialPointsDataFrame-method}
\alias{coerce,SegLocal,SpatialPixelsDataFrame-method}
\alias{coerce,SegLocal,sf-method}
\alias{coerce,SpatialPointsDataFrame,SegLocal-method}
\alias{coerce,SpatialPolygonsDataFrame,SegLocal-method}
\alias{show,SegLocal-method}
\alias{print.SegLocal-method}
\alias{as.data.frame.SegLocal}
\alias{plot.SegLocal-method}
\alias{points.SegLocal-method}
\alias{spplot,SegLocal-method}
\alias{summary.SegLocal-method}
\alias{update.SegLocal-method}

\title{Class SegLocal}
\description{A class to hold the local population composition data.}
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{localenv}}, or the constructor \code{\link{SegLocal}}.
}
\section{Slots}{
  \describe{
    \item{coords}{a numeric matrix or data frame with coordinates (each row is a point).}
    \item{data}{an object of class \code{matrix} containing the population data at each data point. The number of rows in \sQuote{data} should equal the number of points in \sQuote{coords}, and the number of columns should be greater than one (i.e., at least two population groups are required).}
    \item{env}{an object of class \code{matrix} containing the local environment parameters. Must be the same dimensions as \sQuote{data}.}
    \item{geometry}{an optional \code{sfc} geometry object with one feature per row. When present, coercion to \code{sf} returns the local environment values with this geometry.}
    \item{proj4string}{an optional projection string of class \code{CRS}.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "SegLocal", to = "SpatialPoints")}: coerce an object of class \code{SegLocal} to an object of class \code{SpatialPoints}. The points have no attribute data.}
    \item{coerce}{\code{signature(from = "SegLocal", to = "SpatialPointsDataFrame")}: coerce an object of class \code{SegLocal} to an object of class \code{SpatialPointsDataFrame}. The values in the slot \sQuote{env} will be used as the attribute data.}
    \item{coerce}{\code{signature(from = "SegLocal", to = "SpatialPixelsDataFrame")}: coerce an object of class \code{SegLocal} to an object of class \code{SpatialPixelsDataFrame}. The values in the slot \sQuote{env} will be used as the attribute data. May not work when the points are irregularly spaced.}
    \item{coerce}{\code{signature(from = "SegLocal", to = "sf")}: coerce an object of class \code{SegLocal} to an \code{sf} object. The values in the slot \sQuote{env} are used as attributes. Stored geometry is used when available; otherwise point geometry is created from \sQuote{coords}.}
    \item{as.data.frame}{return a data frame containing centroid coordinates \code{x} and \code{y}, followed by the local environment columns. Use \code{columns} to select a subset of local environment columns.}
    \item{coerce}{\code{signature(from = "SpatialPointsDataFrame", to = "SegLocal")}: coerce an object of class \code{SpatialPointsDataFrame} to an object of class \code{SegLocal}.}
    \item{coerce}{\code{signature(from = "SpatialPolygonsDataFrame", to = "SegLocal")}: coerce an object of class \code{SpatialPolygonsDataFrame} to an object of class \code{SegLocal}.}
    \item{show}{\code{signature(object = "SegLocal")}: show the number of points and data columns in an object of class \code{SegLocal}.}
    \item{print}{\code{signature(x = "SegLocal")}: same as \code{show}.}
    \item{plot}{\code{signature(x = "SegLocal")}: draw a plot, or plots, of points in an object of class \code{SegLocal}. Use an optional argument \sQuote{which.col} to specify a column of the data that determines the points' sizes. See the examples below for demonstration.}
    \item{points}{\code{signature(x = "SegLocal")}: draw points in an object of class \code{SegLocal} on an active graphic device.}
    \item{spplot}{\code{signature(obj = "SegLocal")}: coerce an object of class \code{SegLocal} to an object of class \code{SpatialPixelsDataFrame} or \code{SpatialPointsDataFrame} and display it. See \code{help(spplot)} for more details about the graphical parameter arguments.}
    \item{summary}{\code{signature(object = "SegLocal")}: summarise the population compositions of points and local environments in an object of class \code{SegLocal}.}
    \item{update}{\code{signature(object = "SegLocal")}: update an existing object of class \code{SegLocal}.}
  }
}
\author{Seong-Yun Hong}
\seealso{
\code{\link{SegLocal}}, \code{\link{localenv}}
}
\examples{
# creates 100 regularly-spaced data points and 3 population groups
xy <- expand.grid(1:10, 1:10)
colnames(xy) <- c("x", "y")
pop <- matrix(runif(300), ncol = 3)
colnames(pop) <- LETTERS[1:3]

# constructs an object of class 'SegLocal'
v <- SegLocal(coords = as.matrix(xy), data = pop, env = pop)
summary(v)
head(as.data.frame(v))

# updates the map projection information
v <- update(v, proj4string = st_crs("+proj=nzmg +datum=nzgd49"))
summary(v)

# converts the object to sf and plots a local population column
if (requireNamespace("ggplot2", quietly = TRUE)) {
  v_sf <- st_as_sf(v, columns = "A")
  ggplot2::ggplot(v_sf) +
    ggplot2::geom_sf(ggplot2::aes(color = A), size = 2) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_minimal()
}
}
\keyword{classes}

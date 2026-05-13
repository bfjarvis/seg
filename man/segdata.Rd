\name{segdata}
\docType{data}
\alias{segdata}
\title{Hypothetical Patterns of Segregation}
\description{
This data set contains eight different spatial configurations that were used by Morrill (1991) and Wong (1993) to test their segregation measures.
}
\usage{
data(segdata)
}
\format{
An object of class \code{data.frame}. The data set contains 16 columns, representing eight idealised spatial patterns. Each column indicates the following information:

\tabular{rll}{
  [,1] \tab A1 \tab Pattern A, Group 1 \cr
  [,2] \tab A2 \tab Pattern A, Group 2 \cr
  [,3] \tab B1 \tab Pattern B, Group 1 \cr
  [,4] \tab B2 \tab Pattern B, Group 2 \cr
  [,5] \tab C1 \tab Pattern C, Group 1 \cr
  [,6] \tab C2 \tab Pattern C, Group 2 \cr
  [,7] \tab D1 \tab Pattern D, Group 1 \cr
  [,8] \tab D2 \tab Pattern D, Group 2 \cr
  [,9] \tab E1 \tab Pattern E, Group 1 \cr
  [,10] \tab E2 \tab Pattern E, Group 2 \cr
  [,11] \tab F1 \tab Pattern F, Group 1 \cr
  [,12] \tab F2 \tab Pattern F, Group 2 \cr
  [,13] \tab G1 \tab Pattern G, Group 1 \cr
  [,14] \tab G2 \tab Pattern G, Group 2 \cr
  [,15] \tab H1 \tab Pattern H, Group 1 \cr
  [,16] \tab H2 \tab Pattern H, Group 2
}
}
\source{
Morrill, R. L. (1991). On the measure of geographic segregation. \emph{Geography Research Forum}, \bold{11}, 25-36.

Wong, D. W. S. (1993). Spatial Indices of Segregation. \emph{Urban Studies}, \bold{30}, 559-572. 
}
\examples{
data(segdata)
geom <- st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)))) |>
  st_make_grid(cellsize = 1)
grd_sf <- st_sf(segdata, geometry = geom)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_data <- do.call(
    rbind,
    lapply(1:8, function(i) {
      idx <- 2 * i - 1
      data.frame(
        pattern = LETTERS[i],
        group_share = segdata[, idx],
        geometry = st_geometry(grd_sf)
      )
    })
  )
  plot_data <- st_sf(plot_data)

  ggplot2::ggplot(plot_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = group_share), color = "white") +
    ggplot2::facet_wrap(~ pattern) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_void()
}
}
\keyword{datasets}

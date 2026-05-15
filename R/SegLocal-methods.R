# ------------------------------------------------------------------------------
# Methods for class 'SegLocal'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Coercion methods
# ------------------------------------------------------------------------------
.seg_env_subset <- function(env, columns = NULL) {
  if (is.null(columns))
    return(env)
  if (is.numeric(columns))
    return(env[, columns, drop = FALSE])
  if (!all(columns %in% colnames(env)))
    stop("'columns' must name columns in the local environment matrix",
         call. = FALSE)
  env[, columns, drop = FALSE]
}

.seg_env_sf <- function(env, coords, geometry, crs, columns = NULL) {
  out <- data.frame(.seg_env_subset(env, columns))
  if (!is.null(geometry)) {
    st_sf(out, geometry = geometry)
  } else {
    out$x <- coords[, 1]
    out$y <- coords[, 2]
    st_as_sf(out, coords = c("x", "y"), crs = st_crs(crs))
  }
}

st_as_sf.SegLocal <- function(x, ..., columns = NULL) {
  validObject(x)
  .seg_env_sf(x@env, x@coords, x@geometry, x@proj4string, columns = columns)
}

as.data.frame.SegLocal <- function(x, row.names = NULL, optional = FALSE, ...,
                                   columns = NULL) {
  validObject(x)
  out <- data.frame(x = x@coords[, 1], y = x@coords[, 2],
                    .seg_env_subset(x@env, columns),
                    check.names = !optional)
  if (!is.null(row.names))
    rownames(out) <- row.names
  out
}

setAs("SegLocal", "sf", 
      function(from) {
        st_as_sf(from)
      })

setAs("SegLocal", "SpatialPointsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPointsDataFrame(coords = from@coords,
                               data = data.frame(from@env),
                               proj4string = as(from@proj4string, "CRS"))
      })

setAs("SegLocal", "SpatialPixelsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPixelsDataFrame(points = from@coords,
                               data = data.frame(from@env),
                               proj4string = as(from@proj4string, "CRS"))
      })

setAs("SpatialPointsDataFrame", "SegLocal", 
      function(from) {
        SegLocal(coords = st_coordinates(from), data = from@data, env = from@data, 
                 proj4string = st_crs(from), geometry = st_geometry(st_as_sf(from)))
      })

setAs("SpatialPolygonsDataFrame", "SegLocal", 
      function(from) {
        from_sf <- st_as_sf(from)
        coords <- st_geometry(from_sf) |>
          st_make_valid() |>
          st_centroid() |>
          st_coordinates()
        SegLocal(coords = coords, data = from@data, env = from@data,
                 proj4string = st_crs(from),
                 geometry = st_geometry(from_sf))
      })

# ------------------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegLocal"), function(object) {
  validObject(object)
  cat("Class                 :", class(object), "\n")
  cat("Number of data points :", nrow(object@coords), "\n")
  cat("Number of data columns:", ncol(object@data), "\n")
  cat("Projection            :", st_crs(object@proj4string)$wkt, "\n")
  cat("Slot names            :", slotNames(object), "\n")
})

print.SegLocal <- function(x, ...) {
  validObject(x)
  cat("Class                 :", class(x), "\n")
  cat("Number of data points :", nrow(x@coords), "\n")
  cat("Number of data columns:", ncol(x@data), "\n")
  cat("Projection            :", st_crs(x@proj4string)$wkt, "\n")
  cat("Slot names            :", slotNames(x), "\n")
}

# ------------------------------------------------------------------------------
# Plotting (using ggplot2 instead of spplot)
# ------------------------------------------------------------------------------
plot.SegLocal <- function(x, which.col = 1:ncol(x@env), main = NULL, ...) {
  validObject(x)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("package 'ggplot2' is required for plotting", call. = FALSE)
  
  df <- as.data.frame(cbind(x@coords, x@env))
  colnames(df) <- c("x", "y", colnames(x@env))
  
  if (is.null(main)) main <- paste("Data", which.col)
  
  # ggplot(df, aes(x = x, y = y)) +
  #   geom_point(aes(size = df[, which.col]), color = "blue", alpha = 0.6) +
  #   scale_size_continuous(range = c(1, 6)) +
  #   theme_minimal() +
  #   labs(title = main, size = "Value")
}

points.SegLocal <- function(x, which.col = 1, ...) {
  validObject(x)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("package 'ggplot2' is required for plotting", call. = FALSE)
  
  df <- as.data.frame(cbind(x@coords, x@env))
  colnames(df) <- c("x", "y", colnames(x@env))
  
  if (length(which.col) > 1) warning("'which.col' has a length > 1", call. = FALSE)
  
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(size = df[, which.col]), color = "red",
                        alpha = 0.6) +
    ggplot2::scale_size_continuous(range = c(1, 6)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Variable:", colnames(x@env)[which.col]),
                  size = "Value")
}

setMethod("spplot", signature(obj = "SegLocal"), function(obj, ...) {
  validObject(obj)
  spO <- try(as(as(obj, "sf"), "Spatial"), silent = TRUE)
  if (inherits(spO, "try-error"))
    stop("failed to convert 'obj' to a spatial object for spplot",
         call. = FALSE)

  spplot(spO, ...)
})

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
summary.SegLocal <- function(object, ...) {
  validObject(object)
  cat("An object of class \"", class(object), "\"\n", sep = "")
  cat("Coordinates:\n")
  tmp <- apply(object@coords, 2, range)
  print(tmp)
  
  if (is.na(st_crs(object@proj4string)$wkt)) {
    cat("Is projected: FALSE\n")
  } else {
    cat("Is projected: TRUE\n")
    cat("Projection  : ", st_crs(object@proj4string)$wkt, "\n")
  }
  
  cat("\nData values (%):\n")
  tmp <- t(apply(object@data, 1, function(z) z / sum(z))) * 100
  print(apply(tmp, 2, summary, ...))
  
  cat("\nLocal environment composition (%):\n")
  tmp <- t(apply(object@env, 1, function(z) z / sum(z))) * 100
  print(apply(tmp, 2, summary, ...))
}

# ------------------------------------------------------------------------------
# Updating
# ------------------------------------------------------------------------------
update.SegLocal <- function(object, coords, data, env, proj4string, ...) {
  validObject(object)
  dots <- list(...)
  
  if (missing(coords)) coords <- object@coords
  if (missing(data)) data <- object@data
  if (missing(env)) env <- object@env
  if (missing(proj4string)) proj4string <- object@proj4string
  geometry <- if ("geometry" %in% names(dots)) dots$geometry else
    object@geometry
  
  SegLocal(coords, data, env, proj4string, geometry = geometry)
}

# Methods that are not so useful (removed on 23 December 2013) ...
#
# as.list.SegLocal <- function(x, ...) {
#   validObject(x)
#   list(coords = x@coords, data = x@data, env = x@env, 
#        proj4string = x@proj4string)
# }
# 
# setAs("list", "SegLocal",
#       function(from) {
#         if (is.null(from$proj4string))
#           SegLocal(coords = from$coords, data = from$data, env = from$env)
#         else
#           SegLocal(coords = from$coords, data = from$data, env = from$env, 
#                    proj4string = from$proj4string)
#       })
# setAs("SegLocal", "list", 
#       function(from) {
#         validObject(from)
#         list(coords = from@coords, data = from@data, env = from@env, 
#              proj4string = from@proj4string)
#       })
#
# "[[.SegLocal" <- function(i, ...) {
#   validObject(x)
#   slotnames <- slotNames(x)
# 
#   if (is.numeric(i)) {
#     i <- as.integer(i)
#     if (i > length(slotnames))
#       chosen <- NULL
#     else {
#       chosen <- slotnames[i]
#       chosen <- paste("x@", chosen, sep = "")
#       chosen <- eval(parse(text = chosen))
#     }
#   }
#   
#   else if (is.character(i)) {
#     chosen <- paste("x@", i, sep = "")
#     chosen <- eval(parse(text = chosen))
#   }
#   
#   else {
#     chosen <- NULL
#   }
# }

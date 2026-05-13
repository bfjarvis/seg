# ------------------------------------------------------------------------------
# Class 'SegResult'
#
# Unified container for local environments, segregation indices, or both.
# ------------------------------------------------------------------------------
setClass(Class = "SegResult",
         slots = c(coords = "ANY", data = "ANY", env = "ANY",
                   bands = "numeric", indices = "list",
                   measures = "character", weighting = "character",
                   power = "numeric", normalize = "logical",
                   neighbors = "list", proj4string = "crs",
                   call = "ANY"))

setValidity(Class = "SegResult",
            method = function(object) {
              if (!is.null(object@coords) &&
                  (!is.matrix(object@coords) || ncol(object@coords) != 2 ||
                   !is.numeric(object@coords)))
                return("'coords' must be NULL or a numeric matrix with two columns")
              if (!is.null(object@data) &&
                  (!is.matrix(object@data) || ncol(object@data) < 2 ||
                   !is.numeric(object@data)))
                return("'data' must be NULL or a numeric matrix with at least two columns")
              if (!is.null(object@coords) && !is.null(object@data) &&
                  nrow(object@coords) != nrow(object@data))
                return("'data' must have the same number of rows as 'coords'")
              if (!is.null(object@env)) {
                if (!is.list(object@env))
                  return("'env' must be NULL or a list of matrices")
                if (length(object@env) != length(object@bands))
                  return("'env' must have one element per band")
              }
              if (length(object@bands) == 0 || any(!is.finite(object@bands)))
                return("'bands' must contain finite numeric values")
              if (class(object@proj4string) != "crs")
                return("'proj4string' is not a valid CRS object")
              TRUE
            })

SegResult <- function(coords = NULL, data = NULL, env = NULL, bands,
                      indices = list(), measures = character(),
                      weighting = character(), power = numeric(),
                      normalize = logical(), neighbors = list(type = "radius"),
                      proj4string = st_crs(as.character(NA)), call = NULL) {
  new("SegResult", coords = coords, data = data, env = env, bands = bands,
      indices = indices, measures = measures, weighting = weighting,
      power = power, normalize = normalize, neighbors = neighbors,
      proj4string = proj4string, call = call)
}

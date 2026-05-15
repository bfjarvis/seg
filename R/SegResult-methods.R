# ------------------------------------------------------------------------------
# Methods for class 'SegResult'
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegResult"), function(object) {
  validObject(object)
  cat("\n\tSpatial segregation result\n\n")
  cat("Bands                 :", length(object@bands), "\n")
  cat("Neighborhood type     :", object@neighbors$type %||% "radius", "\n")
  cat("Local environments    :", if (is.null(object@env)) "not stored" else "stored", "\n")
  cat("Segregation indices   :", if (length(object@indices) == 0) "not stored" else "stored", "\n")
  if (length(object@weighting) > 0) {
    cat("Weighting             :", object@weighting, "\n")
    cat("Power                 :", object@power, "\n")
    cat("Normalize distances   :", object@normalize, "\n")
  }
})

print.SegResult <- function(x, ...) {
  show(x)
}

as.list.SegResult <- function(x, ...) {
  validObject(x)
  list(bands = x@bands, indices = x@indices, env = x@env,
       coords = x@coords, data = x@data)
}

as.data.frame.SegResult <- function(x, row.names = NULL, optional = FALSE, ...,
                                    what = c("indices", "pairwise", "exposure")) {
  validObject(x)
  what <- match.arg(what)
  if (length(x@indices) == 0)
    stop("'x' does not contain segregation indices", call. = FALSE)

  if (what == "indices") {
    overall <- x@indices$overall
    out <- data.frame()
    scalar_names <- intersect(c("d", "r", "h"), names(overall))
    for (nm in scalar_names) {
      values <- overall[[nm]]
      if (length(values) > 0)
        out <- rbind(out, data.frame(comparison = "overall",
                                     band = x@bands, measure = nm,
                                     value = values))
    }
    rownames(out) <- NULL
    return(out)
  }

  if (what == "pairwise") {
    pairwise <- x@indices$pairwise
    out <- data.frame()
    scalar_names <- intersect(c("d", "r", "h"), names(pairwise))
    for (nm in scalar_names) {
      mats <- pairwise[[nm]]
      if (length(mats) == 0)
        next
      out <- rbind(out, do.call(rbind, lapply(seq_along(mats), function(i) {
        mat <- mats[[i]]
        idx <- which(upper.tri(mat), arr.ind = TRUE)
        data.frame(comparison = "pairwise", band = x@bands[i], measure = nm,
                   group_a = rownames(mat)[idx[, 1]],
                   group_b = colnames(mat)[idx[, 2]],
                   value = mat[idx], row.names = NULL)
      })))
    }
    rownames(out) <- NULL
    return(out)
  }

  p <- x@indices$overall$p
  if (is.null(p) || length(p) == 0)
    return(data.frame())
  out <- do.call(rbind, lapply(seq_along(p), function(i) {
    mat <- p[[i]]
    expand <- expand.grid(from_group = rownames(mat), to_group = colnames(mat),
                          stringsAsFactors = FALSE)
    data.frame(band = x@bands[i], expand, exposure = as.vector(mat))
  }))
  rownames(out) <- NULL
  out
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

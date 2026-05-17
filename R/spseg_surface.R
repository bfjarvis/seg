# ------------------------------------------------------------------------------
# Population surface construction for spseg
# ------------------------------------------------------------------------------

spseg_surface_method <- function(surface) {
  match.arg(surface[1], c("raw", "grid", "pycno"))
}

spseg_surface_config <- function(surface, smoothing, dots) {
  deprecated <- character()

  if (!is.null(smoothing)) {
    deprecated <- c(deprecated, "smoothing")
    old <- if (is.list(smoothing)) smoothing else list(smoothing = smoothing)
    old_method <- old$smoothing %||% "none"
    old_method <- match.arg(old_method, c("none", "raw", "equal", "grid",
                                          "kernel", "pycno"),
                            several.ok = FALSE)
    surface <- switch(old_method,
                      none = "raw",
                      raw = "raw",
                      equal = "grid",
                      grid = "grid",
                      kernel = "pycno",
                      pycno = "pycno")
    old$smoothing <- NULL
    if (!is.null(old$nocol) && is.null(old$ncol))
      old$ncol <- old$nocol
    old$nocol <- NULL
    dots <- c(old, dots)
  }

  if (!is.null(dots$nocol) && is.null(dots$ncol))
    dots$ncol <- dots$nocol
  dots$nocol <- NULL

  if (length(deprecated) > 0)
    warning("Deprecated smoothing argument(s) supplied: ",
            paste(unique(deprecated), collapse = ", "),
            ". Use surface = c(\"raw\", \"grid\", \"pycno\") and pass ",
            "surface-construction options through ....", call. = FALSE)

  list(surface = surface, dots = dots)
}

spseg_surface <- function(x, coords, data, surface, args, verbose) {
  surface <- spseg_surface_method(surface)
  if (surface == "raw")
    return(list(coords = coords, data = data, geometry = NULL))
  if (!inherits(x, "sf"))
    stop("surface = \"", surface, "\" requires 'x' to be an sf object",
         call. = FALSE)
  if (surface == "grid")
    return(spseg_surface_grid(x, data, args, verbose))
  spseg_surface_pycno(x, data, args, verbose)
}

spseg_surface_grid <- function(x, data, args, verbose) {
  grid_dims <- spseg_surface_grid_dims(x, args)

  if (verbose) {
    begTime <- Sys.time()
    message("spseg_surface_grid: redistributing population data ...")
  }

  bbox <- st_bbox(x)
  xx <- seq(bbox$xmin, bbox$xmax, length.out = grid_dims$ncol)
  yy <- seq(bbox$ymin, bbox$ymax, length.out = grid_dims$nrow)
  coords <- expand.grid(x = xx, y = yy)
  cell_width <- if (length(xx) > 1) diff(xx)[1] else spseg_surface_celldim(x, args)
  cell_height <- if (length(yy) > 1) diff(yy)[1] else cell_width

  points <- st_as_sf(coords, coords = c("x", "y"), crs = st_crs(x))
  polygon_ids <- st_intersects(points, x, sparse = FALSE)
  polygon_ids <- apply(polygon_ids, 1, function(z) {
    if (any(z)) which(z)[1] else NA_integer_
  })

  keep <- !is.na(polygon_ids)
  coords <- coords[keep, , drop = FALSE]
  polygon_ids <- polygon_ids[keep]
  spseg_surface_warn_missing_zones(polygon_ids, data, "grid")
  geometry <- spseg_surface_cell_geometry(coords, st_crs(x), cell_width,
                                          cell_height)

  cell_per_polygon <- table(polygon_ids)
  values <- matrix(NA_real_, nrow = nrow(coords), ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    counts <- data[polygon_ids, i]
    values[, i] <- counts / cell_per_polygon[as.character(polygon_ids)]
  }

  if (verbose) {
    tt <- as.numeric(difftime(Sys.time(), begTime, units = "sec"))
    message("spseg_surface_grid: done! [", tt, " seconds]")
  }

  colnames(coords) <- c("x", "y")
  colnames(values) <- colnames(data)
  list(coords = coords, data = values, geometry = geometry, id = polygon_ids)
}

spseg_surface_pycno <- function(x, data, args, verbose) {
  grid_dims <- spseg_surface_grid_dims(x, args)
  max_iter <- args$max_iter %||% 1000
  converge <- args$converge %||% 3

  if (verbose) {
    begTime <- Sys.time()
    message("spseg_surface_pycno: pycnophylactic interpolation ...")
  }

  grid <- spseg_surface_grid_index(x, grid_dims)
  coords <- grid$coords
  zone_ids <- grid$zone_ids
  spseg_surface_warn_missing_zones(zone_ids, data, "pycno")
  values <- seg_pycno_cpp(
    zone_ids = zone_ids,
    pops = data,
    nx = grid_dims$ncol,
    ny = grid_dims$nrow,
    max_iter = max_iter,
    converge = converge
  )

  keep <- !is.na(zone_ids)
  coords <- coords[keep, , drop = FALSE]
  values <- values[keep, , drop = FALSE]
  values[is.na(values)] <- 0
  geometry <- spseg_surface_cell_geometry(coords, st_crs(x),
                                          grid$cell_width,
                                          grid$cell_height)
  colnames(coords) <- c("x", "y")
  colnames(values) <- colnames(data)

  if (verbose) {
    tt <- as.numeric(difftime(Sys.time(), begTime, units = "sec"))
    message("spseg_surface_pycno: done! [", tt, " seconds]")
  }

  list(coords = coords, data = values, geometry = geometry)
}

spseg_surface_polygon_ids <- function(coords, x) {
  points <- st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]),
                     coords = c("x", "y"), crs = st_crs(x))
  hits <- st_intersects(points, x, sparse = FALSE)
  apply(hits, 1, function(z) {
    if (any(z)) which(z)[1] else NA_integer_
  })
}

spseg_surface_warn_missing_zones <- function(zone_ids, data, surface) {
  represented <- unique(zone_ids[!is.na(zone_ids)])
  positive <- which(rowSums(data) > 0)
  missing <- setdiff(positive, represented)
  if (length(missing) > 0) {
    warning(
      "surface = \"", surface, "\" did not create grid cells for ",
      length(missing), " non-empty source polygon(s). ",
      "The grid is too coarse to preserve the input population surface; ",
      "use a smaller 'cellsize'/'celldim'.",
      call. = FALSE
    )
  }
}

spseg_surface_grid_index <- function(x, grid_dims) {
  bbox <- st_bbox(x)
  xx <- seq(bbox$xmin, bbox$xmax, length.out = grid_dims$ncol)
  yy <- seq(bbox$ymin, bbox$ymax, length.out = grid_dims$nrow)
  coords <- expand.grid(x = xx, y = yy)
  list(
    coords = coords,
    zone_ids = spseg_surface_polygon_ids(coords, x),
    cell_width = if (length(xx) > 1) diff(xx)[1] else bbox$xmax - bbox$xmin,
    cell_height = if (length(yy) > 1) diff(yy)[1] else bbox$ymax - bbox$ymin
  )
}

spseg_surface_grid_dims <- function(x, args) {
  if (!is.null(args$nrow) && !is.null(args$ncol))
    return(list(nrow = args$nrow, ncol = args$ncol))

  celldim <- spseg_surface_celldim(x, args)
  bbox <- st_bbox(x)
  list(
    nrow = max(1, ceiling((bbox$ymax - bbox$ymin) / celldim) + 1),
    ncol = max(1, ceiling((bbox$xmax - bbox$xmin) / celldim) + 1)
  )
}

spseg_surface_celldim <- function(x, args) {
  celldim <- args$celldim %||% args$cellsize
  if (!is.null(celldim)) {
    if (!is.numeric(celldim) || length(celldim) != 1 || celldim <= 0)
      stop("'celldim' must be a positive numeric scalar", call. = FALSE)
    return(celldim)
  }

  bbox <- st_bbox(x)
  min(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 100
}

spseg_surface_cell_geometry <- function(coords, crs, width, height = width) {
  if (nrow(coords) == 0)
    return(st_sfc(crs = crs))

  cells <- lapply(seq_len(nrow(coords)), function(i) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    st_polygon(list(rbind(
      c(x - width / 2, y - height / 2),
      c(x + width / 2, y - height / 2),
      c(x + width / 2, y + height / 2),
      c(x - width / 2, y + height / 2),
      c(x - width / 2, y - height / 2)
    )))
  })
  st_sfc(cells, crs = crs)
}

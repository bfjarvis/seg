# ------------------------------------------------------------------------------
# Population surface construction for spseg
# ------------------------------------------------------------------------------

#' Construct a population surface for spatial segregation calculations
#'
#' `spseg_surface()` prepares the coordinate and count matrix that will be passed
#' to the segregation engine. The `"raw"` surface returns the input coordinates
#' and data unchanged. The `"grid"` and `"pycno"` surfaces require polygonal
#' `sf` input and redistribute polygon counts to a regular square grid before
#' optionally applying pycnophylactic smoothing.
#'
#' @param x Input object supplied to `spseg()`. Grid-based surfaces require an
#'   `sf` object with polygon geometry.
#' @param coords Numeric two-column coordinate matrix extracted from `x`.
#' @param data Numeric matrix of group counts.
#' @param surface Character scalar. One of `"raw"`, `"grid"`, or `"pycno"`.
#' @param args Named list of surface options passed through `...` from `spseg()`.
#' @param audit Logical. If `TRUE`, retain surface-construction intermediates.
#' @param verbose Logical. If `TRUE`, report timing messages.
#'
#' @return A list with at least `coords`, `data`, and `geometry`. Grid-based
#'   surfaces may also include source-polygon identifiers and fallback metadata.
#'
#' @noRd
spseg_surface <- function(
  x,
  coords,
  data,
  surface,
  args,
  audit = FALSE,
  verbose = FALSE
) {
  if (surface == "raw") {
    return(list(
      coords = coords,
      data = data,
      geometry = NULL,
      surface_info = list(
        type = "raw",
        fallback = list(count = 0L, ids = integer(), cells = integer()),
        audit = NULL
      )
    ))
  }
  if (!inherits(x, "sf")) {
    stop(
      "surface = \"",
      surface,
      "\" requires 'x' to be an sf object",
      call. = FALSE
    )
  }
  if (surface == "grid") {
    return(spseg_surface_grid(x, data, args, audit, verbose))
  }
  spseg_surface_pycno(x, data, args, audit, verbose)
}

#' Redistribute polygon counts to a square grid
#'
#' Builds the gridded population surface used by `surface = "grid"`. Counts are
#' assigned to cells whose centroids fall within each source polygon, with a
#' fallback assignment for non-empty polygons too small to contain a grid-cell
#' centroid.
#'
#' @inheritParams spseg_surface
#'
#' @return A compact surface list with grid-cell coordinates, redistributed group
#'   counts, optional grid geometry, source polygon identifiers, and fallback
#'   metadata.
#'
#' @noRd
spseg_surface_grid <- function(x, data, args, audit = FALSE, verbose = FALSE) {
  if (verbose) {
    beg_time <- Sys.time()
    message("spseg_surface_grid: redistributing population data ...")
  }

  grid <- surface_grid(x, data, args, audit)
  surface_warn_fallback(grid$fallback, "grid")
  out <- surface_compact_grid(
    grid,
    st_crs(x),
    surface_geometry_type(args),
    "grid"
  )

  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), beg_time, units = "sec"))
    message("spseg_surface_grid: done! [", elapsed, " seconds]")
  }

  out
}

#' Build a pycnophylactic population surface
#'
#' First constructs the same count-preserving grid as `spseg_surface_grid()`,
#' then applies the package's internal pycnophylactic smoother. The smoother
#' iteratively averages neighboring cells and corrects each source zone so that
#' group counts remain equal to the gridded counts supplied to the smoother.
#'
#' @inheritParams spseg_surface
#'
#' @return A compact surface list with smoothed grid-cell coordinates, group
#'   counts, optional grid geometry, source polygon identifiers, and fallback
#'   metadata.
#'
#' @noRd
spseg_surface_pycno <- function(x, data, args, audit = FALSE, verbose = FALSE) {
  max_iter <- args$max_iter %||% 1000
  converge <- args$converge %||% 3

  if (verbose) {
    beg_time <- Sys.time()
    message("spseg_surface_pycno: pycnophylactic interpolation ...")
  }

  grid <- surface_grid(x, data, args, audit)
  surface_warn_fallback(grid$fallback, "pycno")
  group_names <- colnames(grid$values)
  grid$values <- seg_pycno_cpp(
    zone_ids = grid$zone_ids,
    values = grid$values,
    nx = grid$ncol,
    ny = grid$nrow,
    max_iter = max_iter,
    converge = converge
  )
  grid$values[is.na(grid$values)] <- 0
  colnames(grid$values) <- group_names
  out <- surface_compact_grid(
    grid,
    st_crs(x),
    surface_geometry_type(args),
    "pycno"
  )

  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), beg_time, units = "sec"))
    message("spseg_surface_pycno: done! [", elapsed, " seconds]")
  }

  out
}

#' Build the full raster-backed grid representation
#'
#' Creates a terra raster template, rasterizes source polygon identifiers,
#' distributes polygon counts over assigned grid cells, and adds the counts of
#' non-empty polygons without a centroid-based grid cell to the cell containing
#' a representative interior point.
#'
#' @param x Polygonal `sf` object.
#' @param data Numeric matrix of group counts, one row per source polygon.
#' @param args Named list of surface options.
#' @param audit Logical. If `TRUE`, retain source and initial-grid information.
#'
#' @return A list describing the full grid, including all cell coordinates,
#'   count values, source-zone identifiers, fallback assignments, grid dimensions,
#'   and cell size.
#'
#' @noRd
surface_grid <- function(x, data, args, audit = FALSE) {
  grid_dims <- surface_grid_dims(x, args)
  template <- surface_raster_template(x, grid_dims)
  zone_ids <- surface_rasterize_zones(x, template)
  fallback <- surface_fallback_counts(zone_ids, x, data, template)
  values <- surface_distribute_counts(zone_ids, data)

  if (length(fallback$cells) > 0) {
    fallback_values <- rowsum(fallback$data, fallback$cells, reorder = FALSE)
    fallback_cells <- as.integer(rownames(fallback_values))
    values[fallback_cells, ] <- values[fallback_cells, , drop = FALSE] +
      fallback_values

    unassigned_cells <- fallback_cells[is.na(zone_ids[fallback_cells])]
    for (cell in unassigned_cells) {
      zone_ids[cell] <- fallback$ids[match(cell, fallback$cells)]
    }
  }

  surface_warn_unrepresented_zones(zone_ids, data, fallback$ids)

  out <- list(
    coords = terra::xyFromCell(template, seq_len(terra::ncell(template))),
    values = values,
    zone_ids = zone_ids,
    fallback = fallback,
    nrow = terra::nrow(template),
    ncol = terra::ncol(template),
    cell_width = terra::xres(template),
    cell_height = terra::yres(template)
  )
  if (audit) {
    out$initial_values <- values
    out$source_data <- data
    out$source_geometry <- st_geometry(x)
  }
  out
}

#' Determine grid dimensions for a surface
#'
#' Uses explicit `nrow` and `ncol` values when supplied. Otherwise, derives grid
#' dimensions from `celldim`/`cellsize` and the input bounding box.
#'
#' @param x Polygonal `sf` object.
#' @param args Named list of surface options.
#'
#' @return A list with integer `nrow` and `ncol` entries.
#'
#' @noRd
surface_grid_dims <- function(x, args) {
  if (!is.null(args$nrow) && !is.null(args$ncol)) {
    return(list(nrow = args$nrow, ncol = args$ncol))
  }
  
  bbox <- st_bbox(x)
  celldim <- args$celldim %||% args$cellsize
  if (is.null(celldim)) {
    celldim <- min(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 100
  }
  if (!is.numeric(celldim) || length(celldim) != 1 || celldim <= 0) {
    stop("'celldim' must be a positive numeric scalar", call. = FALSE)
  }
  
  list(
    nrow = max(1, ceiling((bbox$ymax - bbox$ymin) / celldim)),
    ncol = max(1, ceiling((bbox$xmax - bbox$xmin) / celldim))
  )
}

#' Create a terra raster template for a surface grid
#'
#' @param x Polygonal `sf` object.
#' @param grid_dims List with integer `nrow` and `ncol` entries.
#'
#' @return A terra raster with the extent and CRS of `x`.
#'
#' @noRd
surface_raster_template <- function(x, grid_dims) {
  bbox <- st_bbox(x)
  terra::rast(
    nrows = grid_dims$nrow,
    ncols = grid_dims$ncol,
    xmin = bbox$xmin,
    xmax = bbox$xmax,
    ymin = bbox$ymin,
    ymax = bbox$ymax,
    crs = st_crs(x)$wkt
  )
}

#' Rasterize source polygon identifiers
#'
#' Assigns each grid cell to the polygon containing its centroid. The returned
#' vector uses `NA` for cells outside the source polygons.
#'
#' @param x Polygonal `sf` object.
#' @param template Terra raster template.
#'
#' @return Integer vector of source polygon row numbers, one value per raster
#'   cell.
#'
#' @noRd
surface_rasterize_zones <- function(x, template) {
  x$.seg_zone_id <- seq_len(nrow(x))
  zone_raster <- terra::rasterize(
    terra::vect(x),
    template,
    field = ".seg_zone_id",
    background = NA,
    touches = FALSE
  )
  as.integer(terra::values(zone_raster, mat = FALSE))
}

#' Distribute polygon counts over assigned grid cells
#'
#' Divides each source polygon's group counts evenly across all grid cells
#' assigned to that polygon.
#'
#' @param zone_ids Integer vector of source polygon identifiers, one per grid
#'   cell.
#' @param data Numeric matrix of group counts, one row per source polygon.
#'
#' @return Numeric matrix of grid-cell group counts.
#'
#' @noRd
surface_distribute_counts <- function(zone_ids, data) {
  values <- matrix(0, nrow = length(zone_ids), ncol = ncol(data))
  assigned <- !is.na(zone_ids)
  if (!any(assigned)) {
    colnames(values) <- colnames(data)
    return(values)
  }

  ids <- zone_ids[assigned]
  cells_per_polygon <- table(ids)
  for (i in seq_len(ncol(data))) {
    values[assigned, i] <- data[ids, i] / cells_per_polygon[as.character(ids)]
  }
  colnames(values) <- colnames(data)
  values
}

#' Find fallback cells for non-empty polygons without grid-cell centroids
#'
#' Finds source polygons with positive population but no rasterized grid cell,
#' then identifies the grid cell containing a representative interior point.
#' Their counts are subsequently added directly to that grid cell. When several
#' missing polygons select the same cell, their counts are summed.
#'
#' @param zone_ids Integer vector of source polygon identifiers.
#' @param x Polygonal `sf` object.
#' @param data Numeric matrix of group counts.
#' @param template Terra raster template.
#'
#' @return A list with fallback cell numbers, source polygon ids, and count rows.
#'
#' @noRd
surface_fallback_counts <- function(zone_ids, x, data, template) {
  positive <- which(rowSums(data) > 0)
  represented <- unique(zone_ids[!is.na(zone_ids)])
  missing <- setdiff(positive, represented)
  out <- list(
    cells = integer(),
    ids = integer(),
    data = matrix(numeric(), nrow = 0, ncol = ncol(data))
  )
  colnames(out$data) <- colnames(data)
  if (length(missing) == 0) {
    return(out)
  }

  fallback_points <- suppressWarnings(
    st_point_on_surface(st_geometry(x[missing, , drop = FALSE]))
  )
  fallback_xy <- st_coordinates(fallback_points)[, 1:2, drop = FALSE]
  cells <- terra::cellFromXY(template, fallback_xy)
  valid <- !is.na(cells)

  if (any(!valid)) {
    warning(
      "failed to assign fallback grid cells for ",
      sum(!valid),
      " non-empty source polygon(s)",
      call. = FALSE
    )
  }

  out$cells <- as.integer(cells[valid])
  out$ids <- missing[valid]
  out$data <- data[out$ids, , drop = FALSE]
  colnames(out$data) <- colnames(data)
  out
}

#' Drop unused grid cells and attach optional geometry
#'
#' Keeps cells assigned to source polygons and returns the compact
#' representation expected by `spseg()`.
#'
#' @param grid Full grid list from `surface_grid()`.
#' @param crs Coordinate reference system for returned geometry.
#' @param geometry_type Character scalar: `"none"`, `"points"`, or `"polygons"`.
#' @param surface_type Character scalar identifying the constructed surface.
#'
#' @return A compact surface list with coordinates, counts, geometry, source ids,
#'   and fallback metadata.
#'
#' @noRd
surface_compact_grid <- function(grid, crs, geometry_type, surface_type) {
  keep <- !is.na(grid$zone_ids)
  cell_ids <- which(keep)

  coords <- grid$coords[keep, , drop = FALSE]
  values <- grid$values[keep, , drop = FALSE]
  zone_ids <- grid$zone_ids[keep]
  geometry <- surface_make_geometry(
    coords,
    crs,
    geometry_type,
    grid$cell_width,
    grid$cell_height
  )

  colnames(coords) <- c("x", "y")
  colnames(values) <- colnames(grid$values)
  surface_info <- list(
    type = surface_type,
    grid = list(
      nrow = grid$nrow,
      ncol = grid$ncol,
      cell_width = grid$cell_width,
      cell_height = grid$cell_height,
      retained_cells = length(cell_ids)
    ),
    fallback = list(
      count = length(grid$fallback$ids),
      ids = grid$fallback$ids,
      cells = grid$fallback$cells
    ),
    audit = NULL
  )

  if (!is.null(grid$initial_values)) {
    fallback_audit <- data.frame(
      source_id = grid$fallback$ids,
      cell_id = grid$fallback$cells,
      zone_id = grid$zone_ids[grid$fallback$cells]
    )
    if (nrow(fallback_audit) > 0) {
      fallback_audit <- cbind(
        fallback_audit,
        as.data.frame(grid$fallback$data)
      )
    }

    surface_info$audit <- list(
      source = list(
        data = grid$source_data,
        geometry = grid$source_geometry
      ),
      cells = data.frame(
        row = seq_along(cell_ids),
        cell_id = cell_ids,
        zone_id = zone_ids,
        x = coords[, 1],
        y = coords[, 2]
      ),
      fallback = fallback_audit,
      initial_counts = grid$initial_values[keep, , drop = FALSE],
      smoothed_counts = if (identical(surface_type, "pycno")) {
        values
      } else {
        NULL
      }
    )
  }

  list(
    coords = coords,
    data = values,
    geometry = geometry,
    id = zone_ids,
    fallback = list(
      cells = grid$fallback$cells,
      ids = grid$fallback$ids
    ),
    surface_info = surface_info
  )
}

#' Warn when populated source polygons require fallback assignment
#'
#' @param fallback Fallback assignment list from `surface_fallback_counts()`.
#' @param surface Character scalar identifying the surface type.
#'
#' @return Invisibly returns `NULL`; called for its warning side effect.
#'
#' @noRd
surface_warn_fallback <- function(fallback, surface) {
  n <- length(fallback$ids)
  if (n == 0) {
    return(invisible(NULL))
  }

  detail <- if (identical(surface, "pycno")) {
    paste0(
      "These polygons are not preserved as separate zonal constraints during ",
      "pycnophylactic smoothing. "
    )
  } else {
    ""
  }
  warning(
    n,
    " populated source polygon(s) received no grid-cell centroid and ",
    "were assigned to fallback grid cells. ",
    detail,
    "Consider using a finer grid and inspect 'result$surface_info$fallback'.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Warn about non-empty source polygons missing from the grid
#'
#' @param zone_ids Integer vector of source polygon identifiers.
#' @param data Numeric matrix of group counts.
#' @param fallback_ids Integer vector of source polygon ids represented through
#'   fallback assignment.
#'
#' @return Invisibly returns `NULL`; called for its warning side effect.
#'
#' @noRd
surface_warn_unrepresented_zones <- function(
  zone_ids,
  data,
  fallback_ids = integer()
) {
  represented <- unique(c(zone_ids[!is.na(zone_ids)], fallback_ids))
  positive <- which(rowSums(data) > 0)
  missing <- setdiff(positive, represented)
  if (length(missing) > 0) {
    warning(
      "surface construction did not create grid cells for ",
      length(missing),
      " non-empty source polygon(s). ",
      "The grid is too coarse to preserve the input population surface; ",
      "use a smaller 'cellsize'/'celldim'.",
      call. = FALSE
    )
  }
}

#' Normalize the requested surface geometry type
#'
#' Accepts logical and character inputs for `surface_geometry` and maps aliases
#' to the internal geometry labels.
#'
#' @param args Named list of surface options.
#'
#' @return One of `"none"`, `"points"`, or `"polygons"`.
#'
#' @noRd
surface_geometry_type <- function(args) {
  geometry <- args$surface_geometry %||% args$keep_geometry %||% FALSE
  if (is.logical(geometry)) {
    if (length(geometry) != 1 || is.na(geometry)) {
      stop(
        "'surface_geometry' must be TRUE, FALSE, or a geometry type",
        call. = FALSE
      )
    }
    return(if (geometry) "polygons" else "none")
  }

  geometry <- match.arg(
    as.character(geometry),
    c("none", "points", "polygons", "point", "polygon")
  )
  switch(
    geometry,
    point = "points",
    polygon = "polygons",
    geometry
  )
}

#' Build optional grid-cell geometry
#'
#' @param coords Numeric matrix of grid-cell center coordinates.
#' @param crs Coordinate reference system for the returned geometry.
#' @param type Character scalar: `"none"`, `"points"`, or `"polygons"`.
#' @param width Numeric grid-cell width.
#' @param height Numeric grid-cell height.
#'
#' @return An `sfc` geometry vector, or `NULL` when `type = "none"`.
#'
#' @noRd
surface_make_geometry <- function(coords, crs, type, width, height = width) {
  if (identical(type, "none")) {
    return(NULL)
  }
  if (nrow(coords) == 0) {
    return(st_sfc(crs = crs))
  }

  if (identical(type, "points")) {
    points <- as.data.frame(coords)
    names(points) <- c("x", "y")
    return(st_geometry(st_as_sf(points, coords = c("x", "y"), crs = crs)))
  }

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

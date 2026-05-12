# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- function(x, data, method = "all", smoothing = "none", 
                  nrow = 100, ncol = 100, window, sigma, useC = TRUE, negative.rm = FALSE, 
                  verbose = FALSE, ...) {
  dots <- spseg_dots(...)

  if (inherits(x, "SegLocal"))
    return(spseg_from_localenv(x, method, useC, negative.rm))
  
  # ----------------------------------------------------------------------------
  # STEP 1 Data preparation
  # ----------------------------------------------------------------------------
  if (verbose)
    tmp <- chksegdata(x, data)
  else
    tmp <- suppressMessages(chksegdata(x, data))
  
  coords <- tmp$coords
  data <- tmp$data
  proj4string <- tmp$proj4string
  
  tmp <- spseg_surface(x, coords, data, smoothing, nrow, ncol, window, sigma,
                       verbose)
  coords <- tmp$coords
  data <- tmp$data
  
  # ----------------------------------------------------------------------------
  # STEP 3 Calculate the population composition of each local environment
  # ----------------------------------------------------------------------------
  env <- do.call(localenv, c(list(x = coords, data = data), dots))
  env <- update(env, proj4string = st_crs(proj4string))
  
  # ----------------------------------------------------------------------------
  # STEP 4 Compute the segregation indices
  # ----------------------------------------------------------------------------
  spseg_from_localenv(env, method, useC, negative.rm)
}

# ------------------------------------------------------------------------------
# Function 'spatseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spatseg <- function(env, method = "all", useC = TRUE, negative.rm = FALSE, ...) {
  .Deprecated("spseg")
  spseg(env, measures = method, useC = useC, negative.rm = negative.rm,
        output = "legacy", ...)
}

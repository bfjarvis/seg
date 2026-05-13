#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

/* -----------------------------------------------------------------------------
   envconstruct()

   Start    : 19 April 2011
   Updated  : 23 September 2018
   R objects: xval, yval, data, as.integer(dim), power, as.integer(weighting), 
              as.integer(normalize), threshold
----------------------------------------------------------------------------- */
SEXP envconstruct(SEXP x, SEXP y, SEXP v, SEXP dim, SEXP p, SEXP weighting,
                  SEXP normalize, SEXP d)
{
  int i, j, k, nrow, ncol = INTEGER(dim)[0],
      weighting_id = INTEGER(weighting)[0], normalized = INTEGER(normalize)[0],
      *inBand;
  double *xP, *yP, *vP, *envP, weight, weightSum, dx, dy, dxy, bw, local_tol,
         *rowDist,
         dist = REAL(d)[0], power = REAL(p)[0];
  SEXP env;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(y = coerceVector(y, REALSXP));
  PROTECT(v = coerceVector(v, REALSXP));
  
  nrow = length(x);  
  PROTECT(env = allocMatrix(REALSXP, nrow, ncol));

  xP = REAL(x); yP = REAL(y); vP = REAL(v);
  envP = REAL(env);
  rowDist = (double*) R_alloc(nrow, sizeof(double));
  inBand = (int*) R_alloc(nrow, sizeof(int));
  
  // Calculate the local environment values of each point in turn
  for(i = 0; i < nrow; i++) {
    weightSum = 0;
    bw = normalized ? dist : 1.0;
    if (bw <= 0)
      bw = 1.0;
    local_tol = 0.0;

    // Distance calculation for inverse distance case.
    // Use minimum neighbor distance to adjust distance to self.
    if (weighting_id == 2) {
      for (j = 0; j < nrow; j++) {
        inBand[j] = 0;
        rowDist[j] = 0.0;
        dx = fabs(xP[i] - xP[j]); dy = fabs(yP[i] - yP[j]);
        if (dist >= 0) {
          if (dx > dist || dy > dist)
            continue;
          dxy = sqrt(dx*dx + dy*dy);
          if (dxy > dist)
            continue;
        } else {
          dxy = sqrt(dx*dx + dy*dy);
        }
        inBand[j] = 1;
        rowDist[j] = dxy;
        if (dxy > 0 && (local_tol <= 0 || dxy < local_tol))
          local_tol = dxy;
      }
      local_tol = local_tol > 0 ? local_tol / 2.0 : bw;
    }

    // Calculate the euclidean distance between the point "i" and all other
    // points (including itself) i.e., j = {1, 2, 3, ..., nrow}
    for (j = 0; j < nrow; j++) {
      if (weighting_id == 2) {
        if (!inBand[j])
          continue;
        dxy = rowDist[j];
      } else {
        dx = fabs(xP[i] - xP[j]); dy = fabs(yP[i] - yP[j]);
        // Maximum search radius is given. Any points that are further than the
        // specified distance from the current ont (point "i") will be ignored.
        if (dist >= 0) {
          if (dx > dist || dy > dist)
            continue;      
          dxy = sqrt(dx*dx + dy*dy);
          if (dxy > dist)
            continue;
        } 
        // Maximum search radius is not given. Consider all points.
        else {
          dxy = sqrt(dx*dx + dy*dy);
        }
      }
      
      switch(weighting_id) {
        case 0:
          weight = 1.0;
          break;
        case 1:
          weight = pow(1 - pow(dxy / bw, power), power);
          break;
        case 2:
          if (dxy == 0)
            dxy = local_tol;
          weight = 1 / pow((dxy / bw), power);
          break;
        case 3:
          weight = exp(-dxy / bw);
          break;
        default:
          weight = exp(-dxy / bw);
      }

      if (weighting_id == 1 && dxy > bw)
        weight = 0.0;
      if (weight < 0)
        weight = 0.0;

      // Get weighted total
      for (k = 0; k < ncol; k++) {
        if (weightSum == 0)
          envP[i+nrow*k] = weight * vP[j+nrow*k];
        else
          envP[i+nrow*k] += weight * vP[j+nrow*k];
      }
      weightSum += weight;
    }
    
    // Get weighted average
    for (k = 0; k < ncol; k++)
      envP[i+nrow*k] = envP[i+nrow*k] / weightSum;
  }

  UNPROTECT(4);
  return(env);
}

#include <math.h>
#include <R.h>
#include <Rinternals.h>

static double seg_weight(double dxy, double bw, double power, int weighting_id)
{
  double weight;

  switch(weighting_id) {
    case 0:
      weight = 1.0;
      break;
    case 1:
      weight = pow(1 - pow(dxy / bw, power), power);
      if (dxy > bw)
        weight = 0.0;
      break;
    case 2:
      weight = 1 / pow((dxy / bw), power);
      break;
    case 3:
      weight = exp(-dxy / bw);
      break;
    default:
      weight = exp(-dxy / bw);
  }

  return weight < 0 ? 0.0 : weight;
}

SEXP seg_engine(SEXP x, SEXP y, SEXP v, SEXP dim, SEXP bands, SEXP p,
                SEXP weighting, SEXP normalize, SEXP measures,
                SEXP keep_env, SEXP keep_indices)
{
  int i, j, k, b, m, n, nrow, ncol = INTEGER(dim)[0],
      nbands, weighting_id = INTEGER(weighting)[0],
      normalized = INTEGER(normalize)[0], store_env = INTEGER(keep_env)[0],
      store_indices = INTEGER(keep_indices)[0], *measure_flags, INDEX;
  double *xP, *yP, *vP, *bandsP, power = REAL(p)[0], maxband, bw,
         dx, dy, dxy, weight, weightSum, local_tol, envSum, xSum = 0.0,
         logM, E = 0.0, I = 0.0, Ep, Ip, denominator, tmp,
         *rowDist, *envRow, *envProp, *xRowSum, *xColSum, *xProp,
         *dAcc, *rAcc, *hAcc, *pAcc;
  SEXP out, out_names, env_list = R_NilValue, indices, indices_names,
       d_out, r_out, h_out, p_list, env_mat, p_mat;

  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(y = coerceVector(y, REALSXP));
  PROTECT(v = coerceVector(v, REALSXP));
  PROTECT(bands = coerceVector(bands, REALSXP));
  PROTECT(measures = coerceVector(measures, INTSXP));

  xP = REAL(x); yP = REAL(y); vP = REAL(v); bandsP = REAL(bands);
  measure_flags = INTEGER(measures);
  nrow = length(x);
  nbands = length(bands);

  rowDist = (double*) R_alloc(nrow, sizeof(double));
  envRow = (double*) R_alloc(ncol, sizeof(double));
  envProp = (double*) R_alloc(ncol, sizeof(double));
  xRowSum = (double*) R_alloc(nrow, sizeof(double));
  xColSum = (double*) R_alloc(ncol, sizeof(double));
  xProp = (double*) R_alloc(ncol, sizeof(double));
  dAcc = (double*) R_alloc(nbands, sizeof(double));
  rAcc = (double*) R_alloc(nbands, sizeof(double));
  hAcc = (double*) R_alloc(nbands, sizeof(double));
  pAcc = (double*) R_alloc(nbands * ncol * ncol, sizeof(double));

  for (b = 0; b < nbands; b++) {
    dAcc[b] = 0.0;
    rAcc[b] = 0.0;
    hAcc[b] = 0.0;
  }
  for (i = 0; i < nbands * ncol * ncol; i++)
    pAcc[i] = 0.0;

  maxband = bandsP[0];
  for (b = 1; b < nbands; b++)
    if (bandsP[b] > maxband)
      maxband = bandsP[b];

  for (k = 0; k < ncol; k++)
    xColSum[k] = 0.0;

  for (i = 0; i < nrow; i++) {
    tmp = 0.0;
    for (k = 0; k < ncol; k++) {
      INDEX = i + k * nrow;
      tmp += vP[INDEX];
      xColSum[k] += vP[INDEX];
    }
    xRowSum[i] = tmp;
    xSum += tmp;
  }

  for (k = 0; k < ncol; k++)
    xProp[k] = xColSum[k] / xSum;

  logM = log10(ncol);
  for (k = 0; k < ncol; k++) {
    if (xProp[k] > 0)
      E -= xProp[k] * log10(xProp[k]) / logM;
    I += xProp[k] * (1 - xProp[k]);
  }

  if (store_env) {
    PROTECT(env_list = allocVector(VECSXP, nbands));
    for (b = 0; b < nbands; b++) {
      PROTECT(env_mat = allocMatrix(REALSXP, nrow, ncol));
      SET_VECTOR_ELT(env_list, b, env_mat);
      UNPROTECT(1);
    }
  }

  for (i = 0; i < nrow; i++) {
    for (j = 0; j < nrow; j++) {
      rowDist[j] = -1.0;
      dx = fabs(xP[i] - xP[j]);
      dy = fabs(yP[i] - yP[j]);
      if (maxband >= 0 && (dx > maxband || dy > maxband))
        continue;
      dxy = sqrt(dx * dx + dy * dy);
      if (maxband < 0 || dxy <= maxband)
        rowDist[j] = dxy;
    }

    for (b = 0; b < nbands; b++) {
      bw = normalized ? bandsP[b] : 1.0;
      if (bw <= 0)
        bw = 1.0;

      local_tol = 0.0;
      if (weighting_id == 2) {
        for (j = 0; j < nrow; j++) {
          dxy = rowDist[j];
          if (dxy > 0 && (bandsP[b] < 0 || dxy <= bandsP[b]) &&
              (local_tol <= 0 || dxy < local_tol))
            local_tol = dxy;
        }
        local_tol = local_tol > 0 ? local_tol / 2.0 : bw;
      }

      for (k = 0; k < ncol; k++)
        envRow[k] = 0.0;
      weightSum = 0.0;

      for (j = 0; j < nrow; j++) {
        dxy = rowDist[j];
        if (dxy < 0 || (bandsP[b] >= 0 && dxy > bandsP[b]))
          continue;
        if (weighting_id == 2 && dxy == 0)
          dxy = local_tol;

        weight = seg_weight(dxy, bw, power, weighting_id);
        weightSum += weight;
        for (k = 0; k < ncol; k++)
          envRow[k] += weight * vP[j + k * nrow];
      }

      for (k = 0; k < ncol; k++)
        envRow[k] = envRow[k] / weightSum;

      if (store_env) {
        env_mat = VECTOR_ELT(env_list, b);
        for (k = 0; k < ncol; k++)
          REAL(env_mat)[i + k * nrow] = envRow[k];
      }

      if (store_indices) {
        envSum = 0.0;
        for (k = 0; k < ncol; k++)
          envSum += envRow[k];
        for (k = 0; k < ncol; k++)
          envProp[k] = envSum > 0 ? envRow[k] / envSum : 0.0;

        if (measure_flags[0] == 1) {
          for (m = 0; m < ncol; m++) {
            for (n = 0; n < ncol; n++) {
              INDEX = b * ncol * ncol + m * ncol + n;
              pAcc[INDEX] += vP[i + m * nrow] / xColSum[m] * envProp[n];
            }
          }
        }

        if (measure_flags[1] == 1) {
          Ep = 0.0;
          for (k = 0; k < ncol; k++)
            if (envProp[k] > 0)
              Ep -= envProp[k] * log10(envProp[k]) / logM;
          hAcc[b] += xRowSum[i] * Ep;
        }

        if (measure_flags[2] == 1) {
          Ip = 0.0;
          for (k = 0; k < ncol; k++)
            Ip += envProp[k] * (1 - envProp[k]);
          rAcc[b] += xRowSum[i] * Ip;
        }

        if (measure_flags[3] == 1) {
          denominator = 2 * xSum * I;
          for (k = 0; k < ncol; k++)
            dAcc[b] += (xRowSum[i] / denominator) *
              fabs(envProp[k] - xProp[k]);
        }
      }
    }
  }

  PROTECT(indices = allocVector(VECSXP, 4));
  PROTECT(indices_names = allocVector(STRSXP, 4));
  SET_STRING_ELT(indices_names, 0, mkChar("d"));
  SET_STRING_ELT(indices_names, 1, mkChar("r"));
  SET_STRING_ELT(indices_names, 2, mkChar("h"));
  SET_STRING_ELT(indices_names, 3, mkChar("p"));
  setAttrib(indices, R_NamesSymbol, indices_names);

  if (store_indices && measure_flags[3] == 1) {
    PROTECT(d_out = allocVector(REALSXP, nbands));
    for (b = 0; b < nbands; b++)
      REAL(d_out)[b] = dAcc[b];
    SET_VECTOR_ELT(indices, 0, d_out);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(indices, 0, allocVector(REALSXP, 0));
  }

  if (store_indices && measure_flags[2] == 1) {
    PROTECT(r_out = allocVector(REALSXP, nbands));
    for (b = 0; b < nbands; b++)
      REAL(r_out)[b] = 1 - rAcc[b] / (xSum * I);
    SET_VECTOR_ELT(indices, 1, r_out);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(indices, 1, allocVector(REALSXP, 0));
  }

  if (store_indices && measure_flags[1] == 1) {
    PROTECT(h_out = allocVector(REALSXP, nbands));
    for (b = 0; b < nbands; b++)
      REAL(h_out)[b] = 1 - hAcc[b] / (xSum * E);
    SET_VECTOR_ELT(indices, 2, h_out);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(indices, 2, allocVector(REALSXP, 0));
  }

  if (store_indices && measure_flags[0] == 1) {
    PROTECT(p_list = allocVector(VECSXP, nbands));
    for (b = 0; b < nbands; b++) {
      PROTECT(p_mat = allocMatrix(REALSXP, ncol, ncol));
      for (m = 0; m < ncol; m++)
        for (n = 0; n < ncol; n++)
          REAL(p_mat)[m + n * ncol] = pAcc[b * ncol * ncol + m * ncol + n];
      SET_VECTOR_ELT(p_list, b, p_mat);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(indices, 3, p_list);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(indices, 3, allocVector(VECSXP, 0));
  }

  PROTECT(out = allocVector(VECSXP, 2));
  PROTECT(out_names = allocVector(STRSXP, 2));
  SET_STRING_ELT(out_names, 0, mkChar("env"));
  SET_STRING_ELT(out_names, 1, mkChar("indices"));
  setAttrib(out, R_NamesSymbol, out_names);
  SET_VECTOR_ELT(out, 0, store_env ? env_list : R_NilValue);
  SET_VECTOR_ELT(out, 1, indices);

  UNPROTECT(store_env ? 10 : 9);
  return out;
}

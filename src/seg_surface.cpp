#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

arma::vec pycno_targets(const arma::ivec& zone_ids, const arma::vec& values)
{
  int Z = 0;
  for (arma::uword i = 0; i < zone_ids.n_elem; ++i) {
    if (zone_ids[i] >= 0)
      Z = std::max(Z, zone_ids[i] + 1);
  }

  arma::vec targets(Z, arma::fill::zeros);
  for (arma::uword i = 0; i < zone_ids.n_elem; ++i) {
    const int z = zone_ids[i];
    if (z < 0 || !std::isfinite(values[i]))
      continue;
    targets[z] += values[i];
  }

  return targets;
}

arma::vec pycno_correct_zones(const arma::vec& values,
                              const arma::ivec& zone_ids,
                              const arma::vec& pops)
{
  arma::vec out = values;
  const int n = zone_ids.n_elem;
  const int Z = pops.n_elem;
  arma::vec totals(Z, arma::fill::zeros);
  arma::ivec counts(Z, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    const int z = zone_ids[i];
    if (z < 0 || z >= Z)
      continue;
    if (std::isfinite(out[i]))
      totals[z] += out[i];
    counts[z] += 1;
  }

  for (int i = 0; i < n; ++i) {
    const int z = zone_ids[i];
    if (z < 0 || z >= Z) {
      out[i] = NA_REAL;
      continue;
    }

    const double target = pops[z];
    if (target <= 0.0) {
      out[i] = 0.0;
    } else if (std::isfinite(totals[z]) && totals[z] > 0.0) {
      out[i] *= target / totals[z];
    } else {
      out[i] = target / counts[z];
    }
  }

  return out;
}

arma::vec pycno_neighbor_average(const arma::vec& values,
                                 const arma::ivec& zone_ids,
                                 int nx, int ny)
{
  arma::vec out(values.n_elem, arma::fill::value(NA_REAL));

  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x) {
      const int id = x + y * nx;
      if (zone_ids[id] < 0)
        continue;

      double sum = 0.0;
      int count = 0;
      for (int dy = -1; dy <= 1; ++dy) {
        const int yy = y + dy;
        if (yy < 0 || yy >= ny)
          continue;
        for (int dx = -1; dx <= 1; ++dx) {
          const int xx = x + dx;
          if (xx < 0 || xx >= nx)
            continue;
          const int q = xx + yy * nx;
          if (zone_ids[q] < 0)
            continue;
          if (std::isfinite(values[q])) {
            sum += values[q];
            count += 1;
          }
        }
      }
      out[id] = count > 0 ? sum / count : 0.0;
    }
  }

  return out;
}

arma::vec pycno_group(const arma::ivec& zone_ids, const arma::vec& initial,
                      int nx, int ny, int max_iter, double converge)
{
  arma::vec values = initial;
  const arma::vec targets = pycno_targets(zone_ids, initial);
  for (arma::uword i = 0; i < values.n_elem; ++i) {
    if (zone_ids[i] < 0)
      values[i] = NA_REAL;
    else if (!std::isfinite(values[i]))
      values[i] = 0.0;
  }

  double max_value = 0.0;
  for (arma::uword i = 0; i < values.n_elem; ++i) {
    if (zone_ids[i] >= 0 && std::isfinite(values[i]))
      max_value = std::max(max_value, values[i]);
  }
  double stopper = max_value * std::pow(10.0, -converge);
  if (!std::isfinite(stopper) || stopper <= 0.0)
    stopper = std::numeric_limits<double>::epsilon();

  for (int iter = 0; iter < max_iter; ++iter) {
    const arma::vec old = values;
    values = pycno_neighbor_average(values, zone_ids, nx, ny);
    values = pycno_correct_zones(values, zone_ids, targets);

    double change = 0.0;
    for (arma::uword i = 0; i < values.n_elem; ++i) {
      if (zone_ids[i] >= 0 && std::isfinite(old[i]) && std::isfinite(values[i]))
        change = std::max(change, std::abs(old[i] - values[i]));
    }
    if (change < stopper)
      break;
  }

  return values;
}

} // namespace

// [[Rcpp::export]]
NumericMatrix seg_pycno_cpp(IntegerVector zone_ids, NumericMatrix values,
                            int nx, int ny, int max_iter, double converge)
{
  const int n = zone_ids.size();
  const int M = values.ncol();
  if (n != nx * ny)
    stop("'zone_ids' length must equal 'nx * ny'");
  if (values.nrow() != n)
    stop("'values' must have one row per grid cell");

  arma::ivec zones(n);
  for (int i = 0; i < n; ++i)
    zones[i] = zone_ids[i] == NA_INTEGER ? -1 : zone_ids[i] - 1;

  arma::mat V(values.begin(), n, M, false);
  NumericMatrix out(n, M);

  for (int m = 0; m < M; ++m) {
    const arma::vec smoothed = pycno_group(zones, V.col(m), nx, ny, max_iter,
                                           converge);
    for (int i = 0; i < n; ++i)
      out(i, m) = smoothed[i];
  }

  return out;
}

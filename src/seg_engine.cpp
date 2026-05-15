#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

struct GroupPair {
  int a;
  int b;
  double T_ab;
  double P_a;
  double P_b;
  double E_ab;
  double I_ab;
  bool valid;
};

double seg_weight(double d, double bw, double power, int weighting_id)
{
  double w = 0.0;

  switch (weighting_id) {
  case 0:
    w = 1.0;
    break;
  case 1:
    w = std::pow(1.0 - std::pow(d / bw, power), power);
    if (d > bw)
      w = 0.0;
    break;
  case 2:
    w = 1.0 / std::pow(d / bw, power);
    break;
  case 3:
    w = std::exp(-d / bw);
    break;
  default:
    w = std::exp(-d / bw);
  }

  return w < 0.0 ? 0.0 : w;
}

double entropy(const arma::vec& p_m, double log_base)
{
  double out = 0.0;
  for (arma::uword m = 0; m < p_m.n_elem; ++m) {
    if (p_m[m] > 0.0)
      out -= p_m[m] * std::log(p_m[m]) / log_base;
  }
  return out;
}

double interaction(const arma::vec& p_m)
{
  return arma::accu(p_m % (1.0 - p_m));
}

void normalize_composition(const arma::vec& L_m, arma::vec& p_m)
{
  const double L = arma::accu(L_m);
  if (L > 0.0)
    p_m = L_m / L;
  else
    p_m.zeros();
}

std::vector<GroupPair> make_pairs(const arma::vec& T_m)
{
  const int M = T_m.size();
  std::vector<GroupPair> pairs;
  pairs.reserve(M * (M - 1) / 2);
  const double log2 = std::log(2.0);

  for (int a = 0; a < M - 1; ++a) {
    for (int b = a + 1; b < M; ++b) {
      GroupPair pair;
      pair.a = a;
      pair.b = b;
      pair.T_ab = T_m[a] + T_m[b];
      pair.valid = T_m[a] > 0.0 && T_m[b] > 0.0;
      pair.P_a = pair.valid ? T_m[a] / pair.T_ab : NA_REAL;
      pair.P_b = pair.valid ? T_m[b] / pair.T_ab : NA_REAL;
      pair.E_ab = pair.valid ?
        -(pair.P_a * std::log(pair.P_a) + pair.P_b * std::log(pair.P_b)) /
          log2 : NA_REAL;
      pair.I_ab = pair.valid ?
        pair.P_a * (1.0 - pair.P_a) + pair.P_b * (1.0 - pair.P_b) :
        NA_REAL;
      pairs.push_back(pair);
    }
  }

  return pairs;
}

// Nomenclature follows Reardon and O'Sullivan (2004)
// P: exposure/isolation
void accumulate_P_overall(arma::cube& P_acc, int b, int i,
                          const arma::mat& N,
                          const arma::vec& T_m,
                          const arma::vec& p_m)
{
  P_acc.slice(b) += (N.row(i).t() / T_m) * p_m.t();
}

// H: information theory
void accumulate_H_overall(arma::vec& H_acc, int b, double T_i,
                          const arma::vec& p_m, double logM)
{
  H_acc[b] += T_i * entropy(p_m, logM);
}

// R: relative diversity
void accumulate_R_overall(arma::vec& R_acc, int b, double T_i,
                          const arma::vec& p_m)
{
  R_acc[b] += T_i * interaction(p_m);
}

// D: dissimilarity
void accumulate_D_overall(arma::vec& D_acc, int b, double T_i,
                          const arma::vec& p_m,
                          const arma::vec& P_m,
                          double denominator)
{
  D_acc[b] += (T_i / denominator) * arma::accu(arma::abs(p_m - P_m));
}

void accumulate_pairwise(arma::mat& D_pair_acc,
                         arma::mat& R_pair_acc,
                         arma::mat& H_pair_acc,
                         int b, int q, const GroupPair& pair,
                         const arma::mat& N,
                         const arma::vec& L_m, int i,
                         const IntegerVector& measure_flags)
{
  if (!pair.valid)
    return;

  const double T_iab = N(i, pair.a) + N(i, pair.b);
  const double L_ab = L_m[pair.a] + L_m[pair.b];
  if (T_iab <= 0.0 || L_ab <= 0.0)
    return;

  const double p_a = L_m[pair.a] / L_ab;
  const double p_b = L_m[pair.b] / L_ab;
  const double log2 = std::log(2.0);

  if (measure_flags[3] == 1) {
    const double D_den = 2.0 * pair.T_ab * pair.I_ab;
    D_pair_acc(b, q) += (T_iab / D_den) *
      (std::abs(p_a - pair.P_a) + std::abs(p_b - pair.P_b));
  }

  if (measure_flags[2] == 1) {
    const double R_i = p_a * (1.0 - p_a) + p_b * (1.0 - p_b);
    R_pair_acc(b, q) += T_iab * R_i;
  }

  if (measure_flags[1] == 1) {
    double H_i = 0.0;
    if (p_a > 0.0)
      H_i -= p_a * std::log(p_a) / log2;
    if (p_b > 0.0)
      H_i -= p_b * std::log(p_b) / log2;
    H_pair_acc(b, q) += T_iab * H_i;
  }
}

NumericVector finalize_scalar_index(const arma::vec& acc,
                                    double denominator)
{
  NumericVector out(acc.n_elem);
  for (arma::uword i = 0; i < acc.n_elem; ++i)
    out[i] = 1.0 - acc[i] / denominator;
  return out;
}

List finalize_P_overall(const arma::cube& P_acc, int B)
{
  List out(B);
  for (int b = 0; b < B; ++b)
    out[b] = wrap(P_acc.slice(b));
  return out;
}

List finalize_pairwise_matrix(const arma::mat& acc,
                              const std::vector<GroupPair>& pairs, int B,
                              int M, int measure_id)
{
  List out(B);
  for (int b = 0; b < B; ++b) {
    NumericMatrix mat(M, M);
    std::fill(mat.begin(), mat.end(), NA_REAL);
    for (arma::uword q = 0; q < pairs.size(); ++q) {
      const GroupPair& pair = pairs[q];
      double value = NA_REAL;
      if (pair.valid) {
        const double raw = acc(b, q);
        if (measure_id == 1)
          value = 1.0 - raw / (pair.T_ab * pair.E_ab);
        else if (measure_id == 2)
          value = 1.0 - raw / (pair.T_ab * pair.I_ab);
        else if (measure_id == 3)
          value = raw;
      }
      mat(pair.a, pair.b) = value;
      mat(pair.b, pair.a) = value;
    }
    out[b] = mat;
  }
  return out;
}

} // namespace

extern "C" SEXP seg_engine(SEXP x, SEXP y, SEXP data, SEXP bands, SEXP power,
                           SEXP weighting, SEXP normalize, SEXP measures,
                           SEXP comparison, SEXP keep_env,
                           SEXP keep_indices)
{
  BEGIN_RCPP
  // Input coordinates and Armadillo views of them
  NumericVector X_input(x); 
  NumericVector Y_input(y);
  const arma::vec X(X_input.begin(), X_input.size(), false);
  const arma::vec Y(Y_input.begin(), Y_input.size(), false);
  // Input counts (n by M) and Armadillo views
  NumericMatrix N_input(data);         
  arma::mat N(N_input.begin(), N_input.nrow(), N_input.ncol(), false); 
  NumericVector BW(bands);             // bandwidth/search-radius values
  IntegerVector measure_flags(measures); // exposure, information, diversity, dissimilarity
  IntegerVector comparison_flags(comparison); // overall and pairwise switches

  const int n = X.size();              // number of focal units
  const int M = N.n_cols;              // number of groups
  const int B = BW.size();             // number of bandwidths
  const double weight_power = as<double>(power); // kernel power parameter
  const int weighting_id = as<int>(weighting); // kernel identifier from R
  const bool normalized = as<int>(normalize) == 1; // scale distances by bandwidth
  const bool calculate_overall = comparison_flags[0] == 1; // compute overall indices
  const bool calculate_pairwise = comparison_flags[1] == 1; // compute pairwise indices
  const bool store_env = as<int>(keep_env) == 1; // return local environments
  const bool store_indices = as<int>(keep_indices) == 1; // return segregation indices

  arma::vec d_j(n);                    // distances from focal unit i to all j
  arma::vec L_m(M);                    // local environment counts by group
  arma::vec p_m(M);                    // local environment proportions by group
  arma::vec T_i = arma::sum(N, 1);     // total observed population by unit
  arma::vec T_m = arma::trans(arma::sum(N, 0)); // total observed population by group
  const double T = arma::accu(T_i);    // total observed population
  arma::vec P_m = T_m / T;             // regional group proportions
  arma::vec D_acc(B, arma::fill::zeros); // dissimilarity accumulator by band
  arma::vec R_acc(B, arma::fill::zeros); // relative diversity accumulator by band
  arma::vec H_acc(B, arma::fill::zeros); // information theory accumulator by band
  arma::cube P_acc(M, M, B, arma::fill::zeros); // exposure/isolation accumulator by band

  double max_bw = BW[0];
  for (int b = 1; b < B; ++b)
    max_bw = std::max(max_bw, BW[b]);

  const double logM = std::log(static_cast<double>(M));
  const double E = entropy(P_m, logM);
  const double I = interaction(P_m);

  const std::vector<GroupPair> pairs = make_pairs(T_m);
  const int Q = pairs.size();
  arma::mat D_pair_acc(B, Q, arma::fill::zeros);
  arma::mat R_pair_acc(B, Q, arma::fill::zeros);
  arma::mat H_pair_acc(B, Q, arma::fill::zeros);

  List env_list;
  if (store_env) {
    env_list = List(B);
    for (int b = 0; b < B; ++b)
      env_list[b] = NumericMatrix(n, M);
  }

  for (int i = 0; i < n; ++i) {
    d_j.fill(-1.0);
    for (int j = 0; j < n; ++j) {
      const double dx = std::abs(X[i] - X[j]);
      const double dy = std::abs(Y[i] - Y[j]);
      if (max_bw >= 0.0 && (dx > max_bw || dy > max_bw))
        continue;

      const double d = std::sqrt(dx * dx + dy * dy);
      if (max_bw < 0.0 || d <= max_bw)
        d_j[j] = d;
    }

    for (int b = 0; b < B; ++b) {
      double bw = normalized ? BW[b] : 1.0;
      if (bw <= 0.0)
        bw = 1.0;

      double local_tolerance = 0.0;
      if (weighting_id == 2) {
        for (int j = 0; j < n; ++j) {
          const double d = d_j[j];
          const bool in_band = BW[b] < 0.0 || d <= BW[b];
          if (d > 0.0 && in_band &&
              (local_tolerance <= 0.0 || d < local_tolerance)) {
            local_tolerance = d;
          }
        }
        local_tolerance = local_tolerance > 0.0 ? local_tolerance / 2.0 : bw;
      }

      L_m.zeros();
      double W = 0.0;

      for (int j = 0; j < n; ++j) {
        double d = d_j[j];
        if (d < 0.0 || (BW[b] >= 0.0 && d > BW[b]))
          continue;
        if (weighting_id == 2 && d == 0.0)
          d = local_tolerance;

        const double w_j = seg_weight(d, bw, weight_power, weighting_id);
        W += w_j;
        for (int m = 0; m < M; ++m)
          L_m[m] += w_j * N(j, m);
      }

      for (int m = 0; m < M; ++m)
        L_m[m] /= W;

      if (store_env) {
        NumericMatrix env_matrix = env_list[b];
        for (int m = 0; m < M; ++m)
          env_matrix(i, m) = L_m[m];
      }

      if (!store_indices)
        continue;

      normalize_composition(L_m, p_m);

      if (calculate_overall) {
        if (measure_flags[0] == 1)
          accumulate_P_overall(P_acc, b, i, N, T_m, p_m);

        if (measure_flags[1] == 1)
          accumulate_H_overall(H_acc, b, T_i[i], p_m, logM);

        if (measure_flags[2] == 1)
          accumulate_R_overall(R_acc, b, T_i[i], p_m);

        if (measure_flags[3] == 1)
          accumulate_D_overall(D_acc, b, T_i[i], p_m, P_m, 2.0 * T * I);
      }

      if (calculate_pairwise) {
        for (int q = 0; q < Q; ++q) {
          accumulate_pairwise(D_pair_acc, R_pair_acc, H_pair_acc, b, q,
                              pairs[q], N, L_m, i, measure_flags);
        }
      }
    }
  }

  NumericVector D_out;
  NumericVector R_out;
  NumericVector H_out;
  List P_out;

  if (store_indices && calculate_overall && measure_flags[3] == 1) {
    D_out = NumericVector(B);
    for (int b = 0; b < B; ++b)
      D_out[b] = D_acc[b];
  }

  if (store_indices && calculate_overall && measure_flags[2] == 1)
    R_out = finalize_scalar_index(R_acc, T * I);

  if (store_indices && calculate_overall && measure_flags[1] == 1)
    H_out = finalize_scalar_index(H_acc, T * E);

  if (store_indices && calculate_overall && measure_flags[0] == 1)
    P_out = finalize_P_overall(P_acc, B);

  List overall = List::create(
    _["d"] = D_out,
    _["r"] = R_out,
    _["h"] = H_out,
    _["p"] = P_out
  );

  List pairwise = List::create(
    _["d"] = (store_indices && calculate_pairwise && measure_flags[3] == 1) ?
      finalize_pairwise_matrix(D_pair_acc, pairs, B, M, 3) : List(),
    _["r"] = (store_indices && calculate_pairwise && measure_flags[2] == 1) ?
      finalize_pairwise_matrix(R_pair_acc, pairs, B, M, 2) : List(),
    _["h"] = (store_indices && calculate_pairwise && measure_flags[1] == 1) ?
      finalize_pairwise_matrix(H_pair_acc, pairs, B, M, 1) : List()
  );

  List indices = List::create(
    _["overall"] = overall,
    _["pairwise"] = pairwise
  );

  SEXP env_out = store_env ? static_cast<SEXP>(env_list) : R_NilValue;

  return List::create(
    _["env"] = env_out,
    _["indices"] = indices
  );

  END_RCPP
}

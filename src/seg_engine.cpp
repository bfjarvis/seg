#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Vendored nanoflann provenance, license, and checksum are recorded in
// inst/NOTICE.md.
#include "nanoflann.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
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

struct Neighbor {
  int id;
  double d;
};

struct PointCloud {
  const arma::vec& X;
  const arma::vec& Y;

  inline size_t kdtree_get_point_count() const { return X.n_elem; }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    return dim == 0 ? X[idx] : Y[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, PointCloud>,
  PointCloud, 2, size_t>;

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

bool neighbor_distance_less(const Neighbor& a, const Neighbor& b)
{
  if (a.d == b.d)
    return a.id < b.id;
  return a.d < b.d;
}

std::vector<Neighbor> all_neighbors(int i, const arma::vec& X,
                                    const arma::vec& Y)
{
  const int n = X.n_elem;
  std::vector<Neighbor> neighbors;
  neighbors.reserve(n);
  for (int j = 0; j < n; ++j) {
    const double dx = X[i] - X[j];
    const double dy = Y[i] - Y[j];
    neighbors.push_back({j, std::sqrt(dx * dx + dy * dy)});
  }
  std::sort(neighbors.begin(), neighbors.end(), neighbor_distance_less);
  return neighbors;
}

std::vector<Neighbor> radius_neighbors(const KDTree& tree, int i,
                                       const arma::vec& X,
                                       const arma::vec& Y,
                                       double radius)
{
  if (radius < 0.0)
    return all_neighbors(i, X, Y);

  double query[2] = {X[i], Y[i]};
  std::vector<nanoflann::ResultItem<size_t, double>> results;
  const double radius2 = std::nextafter(radius * radius,
                                        std::numeric_limits<double>::infinity());
  tree.radiusSearch(query, radius2, results,
                    nanoflann::SearchParameters(0.0, true));

  std::vector<Neighbor> neighbors;
  neighbors.reserve(results.size());
  for (const auto& result : results)
    neighbors.push_back({static_cast<int>(result.first),
                         std::sqrt(result.second)});

  std::sort(neighbors.begin(), neighbors.end(), neighbor_distance_less);
  return neighbors;
}

std::vector<Neighbor> knn_neighbors(const KDTree& tree, int i,
                                    const arma::vec& X,
                                    const arma::vec& Y,
                                    const arma::vec& T_i,
                                    double threshold)
{
  const int n = X.n_elem;
  std::vector<Neighbor> neighbors;

  if (T_i[i] >= threshold) {
    neighbors.push_back({i, 0.0});
    return neighbors;
  }

  double query[2] = {X[i], Y[i]};
  int k = std::min(n, 32);

  while (true) {
    std::vector<size_t> ids(k);
    std::vector<double> d2(k);
    const size_t found = tree.knnSearch(query, k, ids.data(), d2.data());

    neighbors.clear();
    neighbors.push_back({i, 0.0});
    for (size_t q = 0; q < found; ++q) {
      const int id = static_cast<int>(ids[q]);
      if (id == i)
        continue;
      neighbors.push_back({id, std::sqrt(d2[q])});
    }
    std::sort(neighbors.begin() + 1, neighbors.end(), neighbor_distance_less);

    double cumulative = T_i[i];
    for (std::size_t q = 1; q < neighbors.size() && cumulative < threshold; ++q)
      cumulative += T_i[neighbors[q].id];

    if (cumulative >= threshold || k >= n)
      break;
    k = std::min(n, k * 2);
  }

  return neighbors;
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
void accumulate_P_multigroup(arma::cube& P_acc, int b, int i,
                          const arma::mat& N,
                          const arma::vec& T_m,
                          const arma::vec& p_m)
{
  P_acc.slice(b) += (N.row(i).t() / T_m) * p_m.t();
}

// H: information theory
void accumulate_H_multigroup(arma::vec& H_acc, int b, double T_i,
                          const arma::vec& p_m, double logM)
{
  H_acc[b] += T_i * entropy(p_m, logM);
}

// R: relative diversity
void accumulate_R_multigroup(arma::vec& R_acc, int b, double T_i,
                          const arma::vec& p_m)
{
  R_acc[b] += T_i * interaction(p_m);
}

// D: dissimilarity
void accumulate_D_multigroup(arma::vec& D_acc, int b, double T_i,
                          const arma::vec& p_m,
                          const arma::vec& P_m,
                          double denominator)
{
  D_acc[b] += (T_i / denominator) * arma::accu(arma::abs(p_m - P_m));
}

void radius_environment(arma::vec& L_m, const std::vector<Neighbor>& neighbors,
                        const arma::mat& N, double band, double power,
                        int weighting_id, bool normalized)
{
  double bw = normalized ? band : 1.0;
  if (bw <= 0.0)
    bw = 1.0;

  double local_tolerance = 0.0;
  if (weighting_id == 2) {
    for (const Neighbor& neighbor : neighbors) {
      const double d = neighbor.d;
      if ((band < 0.0 || d <= band) && d > 0.0 &&
          (local_tolerance <= 0.0 || d < local_tolerance))
        local_tolerance = d;
    }
    local_tolerance = local_tolerance > 0.0 ? local_tolerance / 2.0 : bw;
  }

  L_m.zeros();
  double W = 0.0;
  for (const Neighbor& neighbor : neighbors) {
    double d = neighbor.d;
    if (band >= 0.0 && d > band)
      break;
    if (weighting_id == 2 && d == 0.0)
      d = local_tolerance;

    const double w = seg_weight(d, bw, power, weighting_id);
    W += w;
    L_m += w * N.row(neighbor.id).t();
  }

  if (W > 0.0)
    L_m /= W;
}

void knn_environment(arma::vec& L_m, const std::vector<Neighbor>& neighbors,
                     const arma::mat& N, const arma::vec& T_i,
                     double threshold, double power, int weighting_id,
                     bool normalized)
{
  L_m.zeros();
  if (neighbors.empty())
    return;

  const int focal = neighbors[0].id;
  L_m += N.row(focal).t();
  double cumulative = T_i[focal];
  if (cumulative >= threshold)
    return;

  double farthest = 0.0;
  double local_tolerance = 0.0;
  double remaining = threshold - cumulative;

  for (std::size_t q = 1; q < neighbors.size() && remaining > 0.0; ++q) {
    const int id = neighbors[q].id;
    if (T_i[id] <= 0.0)
      continue;
    const double fraction = std::min(1.0, remaining / T_i[id]);
    if (fraction <= 0.0)
      break;
    farthest = std::max(farthest, neighbors[q].d);
    if (neighbors[q].d > 0.0 &&
        (local_tolerance <= 0.0 || neighbors[q].d < local_tolerance))
      local_tolerance = neighbors[q].d;
    remaining -= fraction * T_i[id];
  }

  double bw = normalized && farthest > 0.0 ?
    farthest * (1.0 + std::sqrt(std::numeric_limits<double>::epsilon())) :
    1.0;
  if (weighting_id == 2)
    local_tolerance = local_tolerance > 0.0 ? local_tolerance / 2.0 : bw;

  remaining = threshold - cumulative;
  for (std::size_t q = 1; q < neighbors.size() && remaining > 0.0; ++q) {
    const int id = neighbors[q].id;
    if (T_i[id] <= 0.0)
      continue;
    const double fraction = std::min(1.0, remaining / T_i[id]);
    if (fraction <= 0.0)
      break;
    double d = neighbors[q].d;
    if (weighting_id == 2 && d == 0.0)
      d = local_tolerance;
    const double w = seg_weight(d, bw, power, weighting_id);
    L_m += fraction * w * N.row(id).t();
    remaining -= fraction * T_i[id];
  }
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

List finalize_P_multigroup(const arma::cube& P_acc, int B)
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

List indices_from_localenvs(const arma::mat& N,
                            const std::vector<arma::mat>& environments,
                            const IntegerVector& measure_flags,
                            const IntegerVector& scope_flags)
{
  const int n = N.n_rows;
  const int M = N.n_cols;
  const int B = environments.size();
  const bool calculate_multigroup = scope_flags[0] == 1;
  const bool calculate_pairwise = scope_flags[1] == 1;

  arma::vec L_m(M);
  arma::vec p_m(M);
  arma::vec T_i = arma::sum(N, 1);
  arma::vec T_m = arma::trans(arma::sum(N, 0));
  const double T = arma::accu(T_i);
  arma::vec P_m = T_m / T;
  arma::vec D_acc(B, arma::fill::zeros);
  arma::vec R_acc(B, arma::fill::zeros);
  arma::vec H_acc(B, arma::fill::zeros);
  arma::cube P_acc(M, M, B, arma::fill::zeros);

  const double logM = std::log(static_cast<double>(M));
  const double E = entropy(P_m, logM);
  const double I = interaction(P_m);

  const std::vector<GroupPair> pairs = make_pairs(T_m);
  const int Q = pairs.size();
  arma::mat D_pair_acc(B, Q, arma::fill::zeros);
  arma::mat R_pair_acc(B, Q, arma::fill::zeros);
  arma::mat H_pair_acc(B, Q, arma::fill::zeros);

  for (int b = 0; b < B; ++b) {
    const arma::mat& L = environments[b];
    for (int i = 0; i < n; ++i) {
      L_m = L.row(i).t();
      normalize_composition(L_m, p_m);

      if (calculate_multigroup) {
        if (measure_flags[0] == 1)
          accumulate_P_multigroup(P_acc, b, i, N, T_m, p_m);

        if (measure_flags[1] == 1)
          accumulate_H_multigroup(H_acc, b, T_i[i], p_m, logM);

        if (measure_flags[2] == 1)
          accumulate_R_multigroup(R_acc, b, T_i[i], p_m);

        if (measure_flags[3] == 1)
          accumulate_D_multigroup(D_acc, b, T_i[i], p_m, P_m, 2.0 * T * I);
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

  if (calculate_multigroup && measure_flags[3] == 1) {
    D_out = NumericVector(B);
    for (int b = 0; b < B; ++b)
      D_out[b] = D_acc[b];
  }

  if (calculate_multigroup && measure_flags[2] == 1)
    R_out = finalize_scalar_index(R_acc, T * I);

  if (calculate_multigroup && measure_flags[1] == 1)
    H_out = finalize_scalar_index(H_acc, T * E);

  if (calculate_multigroup && measure_flags[0] == 1)
    P_out = finalize_P_multigroup(P_acc, B);

  List multigroup = List::create(
    _["d"] = D_out,
    _["r"] = R_out,
    _["h"] = H_out,
    _["p"] = P_out
  );

  List pairwise = List::create(
    _["d"] = (calculate_pairwise && measure_flags[3] == 1) ?
      finalize_pairwise_matrix(D_pair_acc, pairs, B, M, 3) : List(),
    _["r"] = (calculate_pairwise && measure_flags[2] == 1) ?
      finalize_pairwise_matrix(R_pair_acc, pairs, B, M, 2) : List(),
    _["h"] = (calculate_pairwise && measure_flags[1] == 1) ?
      finalize_pairwise_matrix(H_pair_acc, pairs, B, M, 1) : List()
  );

  return List::create(
    _["multigroup"] = multigroup,
    _["pairwise"] = pairwise
  );
}

} // namespace

// [[Rcpp::export]]
List seg_engine_cpp(NumericVector x, NumericVector y, NumericMatrix data,
                    NumericVector bands, double power, int weighting,
                    int normalize, IntegerVector measures,
                    IntegerVector scope, int keep_env, int keep_indices,
                    int neighbors, int search)
{
  // Input coordinates and Armadillo views of them
  const arma::vec X(x.begin(), x.size(), false);
  const arma::vec Y(y.begin(), y.size(), false);
  // Input counts (n by M) and Armadillo views
  arma::mat N(data.begin(), data.nrow(), data.ncol(), false);
  NumericVector BW(bands);             // bandwidth/search-radius values
  IntegerVector measure_flags(measures); // exposure, information, diversity, dissimilarity
  IntegerVector scope_flags(scope); // multigroup and pairwise switches

  const int n = X.size();              // number of focal units
  const int M = N.n_cols;              // number of groups
  const int B = BW.size();             // number of bandwidths
  const double weight_power = power;   // kernel power parameter
  const int weighting_id = weighting;  // kernel identifier from R
  const int neighbors_id = neighbors;  // 0 = radius, 1 = count-based kNN
  const int search_id = search;        // 0 = kd-tree, 1 = brute force
  const bool normalized = normalize == 1; // scale distances by bandwidth
  const bool calculate_multigroup = scope_flags[0] == 1; // compute multigroup indices
  const bool calculate_pairwise = scope_flags[1] == 1; // compute pairwise indices
  const bool store_env = keep_env == 1; // return local environments
  const bool store_indices = keep_indices == 1; // return segregation indices

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

  PointCloud cloud{X, Y};
  std::unique_ptr<KDTree> tree;
  if (search_id == 0)
    tree.reset(new KDTree(
      2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)
    ));

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
    const std::vector<Neighbor> neighbors_i = search_id == 1 ?
      all_neighbors(i, X, Y) :
      (neighbors_id == 1 ?
         knn_neighbors(*tree, i, X, Y, T_i, max_bw) :
         radius_neighbors(*tree, i, X, Y, max_bw));

    for (int b = 0; b < B; ++b) {
      if (neighbors_id == 1) {
        knn_environment(L_m, neighbors_i, N, T_i, BW[b], weight_power,
                        weighting_id, normalized);
      } else {
        radius_environment(L_m, neighbors_i, N, BW[b], weight_power,
                           weighting_id, normalized);
      }

      if (store_env) {
        NumericMatrix env_matrix = env_list[b];
        for (int m = 0; m < M; ++m)
          env_matrix(i, m) = L_m[m];
      }

      if (!store_indices)
        continue;

      normalize_composition(L_m, p_m);

      if (calculate_multigroup) {
        if (measure_flags[0] == 1)
          accumulate_P_multigroup(P_acc, b, i, N, T_m, p_m);

        if (measure_flags[1] == 1)
          accumulate_H_multigroup(H_acc, b, T_i[i], p_m, logM);

        if (measure_flags[2] == 1)
          accumulate_R_multigroup(R_acc, b, T_i[i], p_m);

        if (measure_flags[3] == 1)
          accumulate_D_multigroup(D_acc, b, T_i[i], p_m, P_m, 2.0 * T * I);
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

  if (store_indices && calculate_multigroup && measure_flags[3] == 1) {
    D_out = NumericVector(B);
    for (int b = 0; b < B; ++b)
      D_out[b] = D_acc[b];
  }

  if (store_indices && calculate_multigroup && measure_flags[2] == 1)
    R_out = finalize_scalar_index(R_acc, T * I);

  if (store_indices && calculate_multigroup && measure_flags[1] == 1)
    H_out = finalize_scalar_index(H_acc, T * E);

  if (store_indices && calculate_multigroup && measure_flags[0] == 1)
    P_out = finalize_P_multigroup(P_acc, B);

  List multigroup = List::create(
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
    _["multigroup"] = multigroup,
    _["pairwise"] = pairwise
  );

  SEXP env_out = store_env ? static_cast<SEXP>(env_list) : R_NilValue;

  return List::create(
    _["env"] = env_out,
    _["indices"] = indices
  );
}

// [[Rcpp::export]]
List seg_indices_env_cpp(NumericMatrix data, NumericMatrix env,
                         IntegerVector measures, IntegerVector scope)
{
  arma::mat N(data.begin(), data.nrow(), data.ncol(), false);
  arma::mat L(env.begin(), env.nrow(), env.ncol(), false);
  std::vector<arma::mat> environments;
  environments.push_back(L);

  return indices_from_localenvs(N, environments, measures, scope);
}

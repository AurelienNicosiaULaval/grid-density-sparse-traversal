// sparse_prefix_lbfp.cpp
// Exact sparse-prefix LBFP density evaluation at query points.
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <cmath>
#include <functional>
#include <cstdint>

using namespace Rcpp;
using namespace RcppParallel;

static inline std::string append_key(const std::string& prefix, int value) {
  if (prefix.empty()) return std::to_string(value);
  return prefix + "|" + std::to_string(value);
}

struct PrefixIndex {
  int d;
  std::vector<int> order0;
  std::vector< std::unordered_set<std::string> > prefix_sets;
  std::unordered_map<std::string, int> count_map;

  PrefixIndex(const IntegerMatrix& cells, const IntegerVector& counts, const IntegerVector& order) {
    d = cells.ncol();
    order0.resize(d);
    for (int j = 0; j < d; ++j) order0[j] = order[j] - 1; // R -> C++ indexing
    prefix_sets.resize(d);

    int n_occ = cells.nrow();
    for (int r = 0; r < n_occ; ++r) {
      std::string key = "";
      for (int depth = 0; depth < d; ++depth) {
        int dim = order0[depth];
        key = append_key(key, cells(r, dim));
        prefix_sets[depth].insert(key);
      }
      count_map[key] = counts[r];
    }
  }
};

struct CountIndex {
  int d;
  std::unordered_map<std::string, int> count_map;

  CountIndex(const IntegerMatrix& cells, const IntegerVector& counts) {
    d = cells.ncol();
    int n_occ = cells.nrow();
    for (int r = 0; r < n_occ; ++r) {
      std::string key = "";
      for (int j = 0; j < d; ++j) {
        key = append_key(key, cells(r, j));
      }
      count_map[key] = counts[r];
    }
  }
};

struct LBFPDirectWorker : public Worker {
  const RMatrix<double> X;
  const std::vector<double> origin;
  const std::vector<double> bw;
  const CountIndex& cidx;
  const double denom;
  const std::uint64_t stencil_size;
  RVector<double> density;
  RVector<int> visited;

  LBFPDirectWorker(const NumericMatrix& X_,
                   const std::vector<double>& origin_,
                   const std::vector<double>& bw_,
                   const CountIndex& cidx_,
                   double denom_,
                   std::uint64_t stencil_size_,
                   NumericVector& density_,
                   IntegerVector& visited_)
    : X(X_), origin(origin_), bw(bw_), cidx(cidx_), denom(denom_),
      stencil_size(stencil_size_), density(density_), visited(visited_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int d = cidx.d;
    std::vector<int> idx(d);
    std::vector<double> u(d);

    for (std::size_t i = begin; i < end; ++i) {
      for (int s = 0; s < d; ++s) {
        double z = (X(i, s) - origin[s]) / bw[s];
        int k = static_cast<int>(std::floor(z));
        double us = z - static_cast<double>(k);
        if (us < 0.0 && us > -1e-12) us = 0.0;
        if (us >= 1.0 && us < 1.0 + 1e-12) us = 1.0 - 1e-12;
        if (us < 0.0) us = 0.0;
        if (us >= 1.0) us = 1.0 - 1e-12;
        idx[s] = k;
        u[s] = us;
      }

      double acc = 0.0;
      int vis = 0;

      for (std::uint64_t mask = 0; mask < stencil_size; ++mask) {
        double weight = 1.0;
        std::string key = "";
        for (int s = 0; s < d; ++s) {
          int bit = static_cast<int>((mask >> s) & 1ULL);
          int coord = idx[s] + bit;
          weight *= (bit == 1) ? u[s] : (1.0 - u[s]);
          key = append_key(key, coord);
        }
        auto it = cidx.count_map.find(key);
        if (it != cidx.count_map.end() && it->second != 0) {
          ++vis;
          acc += weight * static_cast<double>(it->second);
        }
      }

      density[i] = acc / denom;
      visited[i] = vis;
    }
  }
};

struct LBFPWorker : public Worker {
  const RMatrix<double> X;
  const std::vector<double> origin;
  const std::vector<double> bw;
  const PrefixIndex& pidx;
  const double denom;
  RVector<double> density;
  RVector<int> visited;
  RVector<int> prefix_nodes;

  LBFPWorker(const NumericMatrix& X_,
             const std::vector<double>& origin_,
             const std::vector<double>& bw_,
             const PrefixIndex& pidx_,
             double denom_,
             NumericVector& density_,
             IntegerVector& visited_,
             IntegerVector& prefix_nodes_)
    : X(X_), origin(origin_), bw(bw_), pidx(pidx_), denom(denom_),
      density(density_), visited(visited_), prefix_nodes(prefix_nodes_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int d = pidx.d;
    std::vector<int> idx(d);
    std::vector<double> u(d);

    for (std::size_t i = begin; i < end; ++i) {
      for (int s = 0; s < d; ++s) {
        double z = (X(i, s) - origin[s]) / bw[s];
        int k = static_cast<int>(std::floor(z));
        double us = z - static_cast<double>(k);
        if (us < 0.0 && us > -1e-12) us = 0.0;
        if (us >= 1.0 && us < 1.0 + 1e-12) us = 1.0 - 1e-12;
        if (us < 0.0) us = 0.0;
        if (us >= 1.0) us = 1.0 - 1e-12;
        idx[s] = k;
        u[s] = us;
      }

      double acc = 0.0;
      int vis = 0;
      int nodes = 0;

      std::function<void(int, const std::string&, double)> rec;
      rec = [&](int depth, const std::string& key, double weight) {
        int dim = pidx.order0[depth];
        for (int bit = 0; bit <= 1; ++bit) {
          int coord = idx[dim] + bit;
          double w = weight * (bit == 1 ? u[dim] : (1.0 - u[dim]));
          std::string new_key = append_key(key, coord);
          ++nodes;
          if (pidx.prefix_sets[depth].find(new_key) == pidx.prefix_sets[depth].end()) continue;
          if (depth == d - 1) {
            auto it = pidx.count_map.find(new_key);
            if (it != pidx.count_map.end() && it->second != 0) {
              ++vis;
              acc += w * static_cast<double>(it->second);
            }
          } else {
            rec(depth + 1, new_key, w);
          }
        }
      };

      rec(0, "", 1.0);
      density[i] = acc / denom;
      visited[i] = vis;
      prefix_nodes[i] = nodes;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List lbfp_sparse_prefix_cpp(Rcpp::NumericMatrix X,
                                  Rcpp::NumericVector origin,
                                  Rcpp::NumericVector bin_width,
                                  int n_train,
                                  Rcpp::IntegerMatrix occupied_cells,
                                  Rcpp::IntegerVector counts,
                                  Rcpp::IntegerVector order) {
  int n_query = X.nrow();
  int d = X.ncol();
  if (n_train <= 0) stop("n_train must be positive.");
  if (origin.size() != d || bin_width.size() != d) stop("origin and bin_width must have length d.");
  if (occupied_cells.ncol() != d) stop("occupied_cells must have d columns.");
  if (counts.size() != occupied_cells.nrow()) stop("counts length must equal occupied_cells rows.");
  if (order.size() != d) stop("order must have length d.");
  std::vector<bool> seen(d, false);
  for (int j = 0; j < d; ++j) {
    if (order[j] < 1 || order[j] > d) stop("order must be a permutation of 1:d.");
    if (seen[order[j] - 1]) stop("order must not contain duplicates.");
    seen[order[j] - 1] = true;
  }

  std::vector<double> origin_v(d), bw_v(d);
  double volume = 1.0;
  for (int j = 0; j < d; ++j) {
    origin_v[j] = origin[j];
    bw_v[j] = bin_width[j];
    if (!R_finite(bw_v[j]) || bw_v[j] <= 0.0) stop("Invalid bin_width.");
    volume *= bw_v[j];
  }
  double denom = static_cast<double>(n_train) * volume;
  PrefixIndex pidx(occupied_cells, counts, order);

  NumericVector density(n_query);
  IntegerVector visited(n_query);
  IntegerVector prefix_nodes(n_query);

  LBFPWorker worker(X, origin_v, bw_v, pidx, denom, density, visited, prefix_nodes);
  parallelFor(0, n_query, worker);

  return List::create(
    _["density"] = density,
    _["visited"] = visited,
    _["prefix_nodes"] = prefix_nodes,
    _["nominal_stencil"] = std::pow(2.0, d)
  );
}

// [[Rcpp::export]]
Rcpp::List lbfp_direct_cpp(Rcpp::NumericMatrix X,
                           Rcpp::NumericVector origin,
                           Rcpp::NumericVector bin_width,
                           int n_train,
                           Rcpp::IntegerMatrix occupied_cells,
                           Rcpp::IntegerVector counts) {
  int n_query = X.nrow();
  int d = X.ncol();
  if (n_train <= 0) stop("n_train must be positive.");
  if (d < 1 || d > 30) stop("Direct LBFP C++ backend requires 1 <= d <= 30.");
  if (origin.size() != d || bin_width.size() != d) stop("origin and bin_width must have length d.");
  if (occupied_cells.ncol() != d) stop("occupied_cells must have d columns.");
  if (counts.size() != occupied_cells.nrow()) stop("counts length must equal occupied_cells rows.");

  std::vector<double> origin_v(d), bw_v(d);
  double volume = 1.0;
  for (int j = 0; j < d; ++j) {
    origin_v[j] = origin[j];
    bw_v[j] = bin_width[j];
    if (!R_finite(bw_v[j]) || bw_v[j] <= 0.0) stop("Invalid bin_width.");
    volume *= bw_v[j];
  }

  std::uint64_t stencil_size = 1ULL << d;
  double denom = static_cast<double>(n_train) * volume;
  CountIndex cidx(occupied_cells, counts);

  NumericVector density(n_query);
  IntegerVector visited(n_query);

  LBFPDirectWorker worker(X, origin_v, bw_v, cidx, denom, stencil_size, density, visited);
  parallelFor(0, n_query, worker);

  return List::create(
    _["density"] = density,
    _["visited"] = visited,
    _["nominal_stencil"] = static_cast<double>(stencil_size)
  );
}

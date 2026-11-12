#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

inline double frac01(double x) {
  double r = std::fmod(x, 1.0);
  if (r < 0.0) r += 1.0;
  return r;
}

double reduced_cost_of_sample(const std::vector<int>& s,
                              const NumericVector& z,
                              const NumericVector& dual,
                              double target) {
  double q = 0.0;
  double yy = 0.0;
  for (int i : s) {
    q += z[i];
    yy += dual[i];
  }
  const double e = q - target;
  return e * e - yy;
}

} // namespace


// [[Rcpp::export]]
List cos_systematic_support_cpp(NumericVector pik,
                                IntegerVector sigma,
                                double tol = 1e-13) {
  const int N = pik.size();
  if (sigma.size() != N) stop("sigma must have length(pik).");

  double sump = 0.0;
  std::vector<unsigned char> seen(N, 0);
  std::vector<int> ord(N);
  for (int k = 0; k < N; ++k) {
    if (!R_finite(pik[k]) || pik[k] <= 0.0 || pik[k] >= 1.0)
      stop("All inclusion probabilities must lie strictly between 0 and 1.");
    sump += pik[k];

    const int u = sigma[k] - 1;
    if (u < 0 || u >= N || seen[u]) stop("sigma must be a permutation of 1,...,N.");
    seen[u] = 1;
    ord[k] = u;
  }

  const int n = static_cast<int>(std::llround(sump));
  if (std::fabs(sump - n) > 1e-8) stop("sum(pik) must be an integer.");

  // Cumulative inclusion probabilities in the supplied order.
  std::vector<double> V(N + 1, 0.0);
  for (int k = 0; k < N; ++k) V[k + 1] = V[k] + pik[ord[k]];

  // Breakpoints on [0,1].  For a generic vector there are N intervals.
  std::vector<double> cuts;
  cuts.reserve(N + 1);
  cuts.push_back(0.0);
  for (int k = 1; k < N; ++k) {
    double r = frac01(V[k]);
    if (r < tol || 1.0 - r < tol) r = 0.0;
    cuts.push_back(r);
  }
  cuts.push_back(1.0);
  std::sort(cuts.begin(), cuts.end());

  std::vector<double> unique_cuts;
  unique_cuts.reserve(cuts.size());
  for (double x : cuts) {
    if (unique_cuts.empty() || x - unique_cuts.back() > tol) {
      unique_cuts.push_back(x);
    } else if (x > unique_cuts.back()) {
      unique_cuts.back() = x;
    }
  }
  if (unique_cuts.front() != 0.0) unique_cuts.insert(unique_cuts.begin(), 0.0);
  if (1.0 - unique_cuts.back() > tol) unique_cuts.push_back(1.0);
  unique_cuts.front() = 0.0;
  unique_cuts.back() = 1.0;

  const int m = static_cast<int>(unique_cuts.size()) - 1;
  IntegerMatrix samples(m, n);
  NumericVector probabilities(m);

  for (int h = 0; h < m; ++h) {
    const double left = unique_cuts[h];
    const double right = unique_cuts[h + 1];
    const double u = 0.5 * (left + right);
    probabilities[h] = right - left;

    int jj = 0;
    for (int k = 0; k < N; ++k) {
      const double a = std::floor(V[k] - u);
      const double b = std::floor(V[k + 1] - u);
      if (b > a) {
        if (jj >= n) stop("Internal error constructing systematic support.");
        samples(h, jj) = ord[k] + 1;  // return population labels, base 1
        ++jj;
      }
    }
    if (jj != n) stop("Internal error: systematic sample has the wrong size.");

    // Canonical representation simplifies duplicate/basis checks in R.
    std::vector<int> tmp(n);
    for (int j = 0; j < n; ++j) tmp[j] = samples(h, j);
    std::sort(tmp.begin(), tmp.end());
    for (int j = 0; j < n; ++j) samples(h, j) = tmp[j];
  }

  return List::create(
    _["samples"] = samples,
    _["probabilities"] = probabilities
  );
}


// [[Rcpp::export]]
List cos_pricing_exact_cpp(NumericVector z,
                           NumericVector dual,
                           double target,
                           int n) {
  const int N = z.size();
  if (dual.size() != N) stop("dual and z must have the same length.");
  if (n <= 0 || n >= N) stop("n must satisfy 0 < n < N.");

  std::vector<int> comb(n);
  for (int j = 0; j < n; ++j) comb[j] = j;

  double best = R_PosInf;
  std::vector<int> best_sample = comb;
  long long evaluated = 0;

  while (true) {
    double rc = reduced_cost_of_sample(comb, z, dual, target);
    ++evaluated;
    if (rc < best) {
      best = rc;
      best_sample = comb;
    }

    int i = n - 1;
    while (i >= 0 && comb[i] == N - n + i) --i;
    if (i < 0) break;

    ++comb[i];
    for (int j = i + 1; j < n; ++j) comb[j] = comb[j - 1] + 1;
  }

  IntegerVector out(n);
  for (int j = 0; j < n; ++j) out[j] = best_sample[j] + 1;

  return List::create(
    _["sample"] = out,
    _["reduced_cost"] = best,
    _["evaluated"] = static_cast<double>(evaluated)
  );
}


// [[Rcpp::export]]
List cos_pricing_swap_cpp(NumericVector z,
                          NumericVector dual,
                          double target,
                          int n,
                          int n_starts = 20,
                          int max_swaps = 2000,
                          double improvement_tol = 1e-12) {
  const int N = z.size();
  if (dual.size() != N) stop("dual and z must have the same length.");
  if (n <= 0 || n >= N) stop("n must satisfy 0 < n < N.");
  if (n_starts < 1) stop("n_starts must be >= 1.");

  RNGScope scope;

  IntegerMatrix local_samples(n_starts, n);
  NumericVector local_values(n_starts);
  IntegerVector local_swaps(n_starts);

  std::vector<int> all(N);
  std::iota(all.begin(), all.end(), 0);

  for (int rr = 0; rr < n_starts; ++rr) {
    // Partial Fisher-Yates shuffle: first n entries form a SRSWOR start.
    std::iota(all.begin(), all.end(), 0);
    for (int i = 0; i < n; ++i) {
      const int j = i + static_cast<int>(std::floor(R::runif(0.0, 1.0) * (N - i)));
      std::swap(all[i], all[j]);
    }

    std::vector<int> selected(n);
    std::vector<unsigned char> in_sample(N, 0);
    double q = 0.0;
    double sy = 0.0;
    for (int i = 0; i < n; ++i) {
      selected[i] = all[i];
      in_sample[selected[i]] = 1;
      q += z[selected[i]];
      sy += dual[selected[i]];
    }

    int accepted = 0;
    while (accepted < max_swaps) {
      const double e = q - target;
      double best_delta = -improvement_tol;
      int best_pos = -1;
      int best_in = -1;

      // Best-improvement exchange neighbourhood.
      for (int a = 0; a < n; ++a) {
        const int i = selected[a];
        const double zi = z[i];
        const double yi = dual[i];

        for (int j = 0; j < N; ++j) {
          if (in_sample[j]) continue;

          const double dz = z[j] - zi;
          const double delta =
            dz * dz + 2.0 * e * dz - dual[j] + yi;

          if (delta < best_delta) {
            best_delta = delta;
            best_pos = a;
            best_in = j;
          }
        }
      }

      if (best_pos < 0) break;

      const int old = selected[best_pos];
      q += z[best_in] - z[old];
      sy += dual[best_in] - dual[old];
      in_sample[old] = 0;
      in_sample[best_in] = 1;
      selected[best_pos] = best_in;
      ++accepted;
    }

    std::sort(selected.begin(), selected.end());
    for (int j = 0; j < n; ++j) local_samples(rr, j) = selected[j] + 1;

    const double ee = q - target;
    local_values[rr] = ee * ee - sy;
    local_swaps[rr] = accepted;
  }

  return List::create(
    _["samples"] = local_samples,
    _["reduced_costs"] = local_values,
    _["accepted_swaps"] = local_swaps
  );
}

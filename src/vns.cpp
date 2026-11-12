#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

inline bool all_equal_pi_std(const std::vector<double>& pi, double tol) {
  if (pi.empty()) return true;
  const double ref = pi[0];
  for (size_t i = 1; i < pi.size(); ++i) {
    if (std::fabs(pi[i] - ref) > tol) return false;
  }
  return true;
}

inline double frac01(double x) {
  double r = std::fmod(x, 1.0);
  if (r < 0.0) r += 1.0;
  return r;
}

inline int gcd_int(int a, int b) {
  a = std::abs(a); b = std::abs(b);
  while (b != 0) { const int r = a % b; a = b; b = r; }
  return a;
}

inline std::vector<double> build_cuts(const std::vector<double>& cumV, double tol) {
  std::vector<double> r;
  r.reserve(cumV.size());
  for (double v : cumV) {
    const double z = frac01(v);
    if (z > tol && z < 1.0 - tol) r.push_back(z);
  }
  std::sort(r.begin(), r.end());
  std::vector<double> cuts;
  cuts.reserve(r.size() + 2);
  cuts.push_back(0.0);
  double last = -1.0;
  for (double z : r) {
    if (last < 0.0 || std::fabs(z - last) > tol) {
      cuts.push_back(z);
      last = z;
    }
  }
  cuts.push_back(1.0);
  return cuts;
}

inline std::vector<int> selected_positions_from_u(const std::vector<double>& cumV, double u) {
  const int N = static_cast<int>(cumV.size());
  std::vector<int> out;
  out.reserve(N);

  // Unit k corresponds to the interval (V_{k-1}, V_k] on the cumulative
  // inclusion-probability line.  It is selected iff that interval contains
  // one of the systematic points u+j.  Testing the change in the integer part
  // is numerically safer than comparing fractional parts.  In particular it
  // treats pi_k = 1 correctly: a certainty unit is selected for every u.
  double prevV = 0.0;
  for (int k = 0; k < N; ++k) {
    const double a = prevV - u;
    const double b = cumV[k] - u;
    if (std::floor(b) > std::floor(a)) out.push_back(k);
    prevV = cumV[k];
  }
  return out;
}

struct BalanceProblem {
  int N;
  int P;
  std::vector<double> pi;
  std::vector<double> inv_pi;
  std::vector<double> X;       // row-major N x P
  std::vector<double> total;   // P
  std::vector<double> scale;   // P, strictly positive
};

BalanceProblem make_problem(const NumericVector& piR,
                            const NumericMatrix& XR,
                            const Nullable<NumericVector>& scaleR) {
  const int N = piR.size();
  const int P = XR.ncol();
  if (XR.nrow() != N) stop("X must have length(pik) rows.");
  if (N < 2) stop("At least two population units are required.");
  if (P < 1) stop("X must have at least one column.");

  BalanceProblem pr;
  pr.N = N;
  pr.P = P;
  pr.pi.assign(piR.begin(), piR.end());
  pr.inv_pi.resize(N);
  pr.X.resize(static_cast<size_t>(N) * P);
  pr.total.assign(P, 0.0);

  long double sump = 0.0L;
  for (int i = 0; i < N; ++i) {
    const double p = pr.pi[i];
    if (!std::isfinite(p) || p <= 0.0 || p > 1.0 + 1e-12) {
      stop("All inclusion probabilities must be in (0, 1].");
    }
    pr.inv_pi[i] = 1.0 / p;
    sump += static_cast<long double>(p);
    for (int j = 0; j < P; ++j) {
      const double x = XR(i, j);
      if (!std::isfinite(x)) stop("X must contain only finite values.");
      pr.X[static_cast<size_t>(i) * P + j] = x;
      pr.total[j] += x;
    }
  }
  if (std::fabs(static_cast<double>(sump - std::round(sump))) > 1e-8) {
    stop("sum(pik) must be an integer for fixed-size systematic sampling.");
  }

  pr.scale.resize(P);
  if (scaleR.isNotNull()) {
    NumericVector sr(scaleR);
    if (sr.size() != P) stop("balance_scale must have ncol(X) entries.");
    for (int j = 0; j < P; ++j) {
      if (!std::isfinite(sr[j]) || sr[j] <= 0.0) stop("balance_scale must be finite and strictly positive.");
      pr.scale[j] = sr[j];
    }
  } else {
    for (int j = 0; j < P; ++j) {
      const double a = std::fabs(pr.total[j]);
      if (a < 1e-14) {
        stop("A column total of X is zero (or numerically zero). Supply balance_scale explicitly.");
      }
      pr.scale[j] = a;
    }
  }
  return pr;
}

std::vector<int> validate_sigma(const IntegerVector& sigmaR, int N) {
  if (sigmaR.size() != N) stop("sigma must have length(pik) entries.");
  std::vector<int> sigma(N);
  std::vector<unsigned char> seen(N, 0);
  for (int k = 0; k < N; ++k) {
    const int u = sigmaR[k] - 1;
    if (u < 0 || u >= N || seen[u]) stop("sigma must be a permutation of 1,...,N.");
    seen[u] = 1;
    sigma[k] = u;
  }
  return sigma;
}

// Exact objective for one systematic design.  Instead of reconstructing the
// sample independently on every support interval, this routine sweeps the
// fractional cumulative-sum breakpoints.  Crossing a breakpoint attached to
// position k moves one systematic point from position k to k+1.  Hence, after
// one O(N) initialization, all support states are updated locally.
struct BreakEvent {
  double r;
  int k;
};

std::vector<BreakEvent> make_break_events(const std::vector<double>& cumV,
                                          double tol_pi) {
  std::vector<BreakEvent> events;
  events.reserve(cumV.size());
  // The last cumulative sum is an integer and corresponds to the circular
  // boundary at zero. It is already absorbed in the initialization.
  for (int k = 0; k + 1 < static_cast<int>(cumV.size()); ++k) {
    const double r = frac01(cumV[k]);
    if (r > tol_pi && r < 1.0 - tol_pi) events.push_back(BreakEvent{r, k});
  }
  std::sort(events.begin(), events.end(), [](const BreakEvent& a, const BreakEvent& b) {
    if (a.r < b.r) return true;
    if (a.r > b.r) return false;
    return a.k < b.k;
  });
  return events;
}

double evaluate_balance_general(const std::vector<int>& sigma,
                                const BalanceProblem& pr,
                                double tol_pi) {
  const int N = pr.N;
  const int P = pr.P;
  std::vector<double> cumV(N);
  if (all_equal_pi_std(pr.pi, tol_pi)) {
    // Avoid cumulative floating-point drift (important when N is large).
    const double p0 = pr.pi[0];
    for (int k = 0; k < N; ++k) cumV[k] = static_cast<double>(k + 1) * p0;
  } else {
    long double acc = 0.0L;
    for (int k = 0; k < N; ++k) {
      acc += static_cast<long double>(pr.pi[sigma[k]]);
      cumV[k] = static_cast<double>(acc);
    }
  }
  const std::vector<BreakEvent> events = make_break_events(cumV, tol_pi);

  // First interval: (0, first breakpoint), or all of (0,1) if there is none.
  const double first_right = events.empty() ? 1.0 : events.front().r;
  const double u0 = 0.5 * first_right;
  const std::vector<int> initial_pos = selected_positions_from_u(cumV, u0);

  std::vector<double> residual(P, 0.0);
  for (int pos : initial_pos) {
    const int unit = sigma[pos];
    const double w = pr.inv_pi[unit];
    const size_t base = static_cast<size_t>(unit) * P;
    for (int j = 0; j < P; ++j) residual[j] += w * pr.X[base + j];
  }
  double score = 0.0;
  for (int j = 0; j < P; ++j) {
    residual[j] = (residual[j] - pr.total[j]) / pr.scale[j];
    score += residual[j] * residual[j];
  }

  double objective = 0.0;
  double left = 0.0;
  size_t idx = 0;
  while (idx < events.size()) {
    const double r = events[idx].r;
    objective += (r - left) * score;

    // Collect all events at the same circular breakpoint.  Consecutive event
    // indices form a chain.  This occurs in particular when one or more
    // intervening units have pi=1: the fractional cumulative sum is unchanged
    // across a certainty unit.  Such certainty units must remain selected on
    // both sides of the breakpoint.  Therefore a chain a,...,b represents one
    // net replacement sigma[a] -> sigma[b+1], not sequential replacements that
    // would incorrectly remove the certainty units from the sample.
    size_t jdx = idx;
    while (jdx < events.size() && std::fabs(events[jdx].r - r) <= tol_pi) ++jdx;

    size_t qevent = idx;
    while (qevent < jdx) {
      int a = events[qevent].k;
      int b = a;
      ++qevent;
      while (qevent < jdx && events[qevent].k == b + 1) {
        b = events[qevent].k;
        ++qevent;
      }

      const int u_out = sigma[a];
      const int u_in = sigma[b + 1];
      const size_t bo = static_cast<size_t>(u_out) * P;
      const size_t bi = static_cast<size_t>(u_in) * P;

      double dot = 0.0;
      double d2 = 0.0;
      for (int qq = 0; qq < P; ++qq) {
        const double d = (pr.inv_pi[u_in] * pr.X[bi + qq] -
                          pr.inv_pi[u_out] * pr.X[bo + qq]) / pr.scale[qq];
        dot += residual[qq] * d;
        d2 += d * d;
      }
      score += 2.0 * dot + d2;
      for (int qq = 0; qq < P; ++qq) {
        const double d = (pr.inv_pi[u_in] * pr.X[bi + qq] -
                          pr.inv_pi[u_out] * pr.X[bo + qq]) / pr.scale[qq];
        residual[qq] += d;
      }
    }
    left = r;
    idx = jdx;
  }
  objective += (1.0 - left) * score;
  return objective;
}

struct EqualTemplate {
  int N;
  int D;
  int B;
  std::vector<double> probs;
  std::vector<std::vector<int>> selected_positions;
  std::vector<uint64_t> position_bits; // N x B bitset of support intervals
};

EqualTemplate make_equal_template(const BalanceProblem& pr, double tol_pi) {
  const int N = pr.N;
  const int n = static_cast<int>(std::llround(pr.pi[0] * static_cast<double>(N)));
  const int g = gcd_int(N, n);
  const int D = N / std::max(g, 1);

  // For equal pik=n/N the distinct breakpoints are exactly h/D.  Constructing
  // them algebraically avoids the severe cumulative-rounding drift that can
  // otherwise create spurious support points for large N.
  std::vector<double> cumV(N);
  const double p0 = static_cast<double>(n) / static_cast<double>(N);
  for (int k = 0; k < N; ++k) cumV[k] = static_cast<double>(k + 1) * p0;

  EqualTemplate tp;
  tp.N = N;
  tp.D = D;
  tp.B = (tp.D + 63) / 64;
  tp.probs.assign(tp.D, 1.0 / static_cast<double>(tp.D));
  tp.selected_positions.resize(tp.D);
  tp.position_bits.assign(static_cast<size_t>(N) * tp.B, uint64_t(0));

  for (int h = 0; h < tp.D; ++h) {
    const double u = (static_cast<double>(h) + 0.5) / static_cast<double>(tp.D);
    tp.selected_positions[h] = selected_positions_from_u(cumV, u);
    const int block = h >> 6;
    const uint64_t bit = uint64_t(1) << (h & 63);
    for (int pos : tp.selected_positions[h]) {
      tp.position_bits[static_cast<size_t>(pos) * tp.B + block] |= bit;
    }
  }
  return tp;
}

struct EqualState {
  std::vector<int> sigma;
  std::vector<double> residual; // D x P, standardized residuals
  std::vector<double> score;    // D
  double f;
};

EqualState make_equal_state(const std::vector<int>& sigma,
                            const BalanceProblem& pr,
                            const EqualTemplate& tp) {
  EqualState st;
  st.sigma = sigma;
  st.residual.assign(static_cast<size_t>(tp.D) * pr.P, 0.0);
  st.score.assign(tp.D, 0.0);
  st.f = 0.0;
  const double w = pr.inv_pi[0];
  std::vector<double> xhat(pr.P);

  for (int h = 0; h < tp.D; ++h) {
    std::fill(xhat.begin(), xhat.end(), 0.0);
    for (int pos : tp.selected_positions[h]) {
      const int unit = st.sigma[pos];
      const size_t base = static_cast<size_t>(unit) * pr.P;
      for (int j = 0; j < pr.P; ++j) xhat[j] += w * pr.X[base + j];
    }
    double sc = 0.0;
    for (int j = 0; j < pr.P; ++j) {
      const double r = (xhat[j] - pr.total[j]) / pr.scale[j];
      st.residual[static_cast<size_t>(h) * pr.P + j] = r;
      sc += r * r;
    }
    st.score[h] = sc;
    st.f += tp.probs[h] * sc;
  }
  return st;
}

inline bool bit_is_set(const EqualTemplate& tp, int pos, int h) {
  const int block = h >> 6;
  const uint64_t bit = uint64_t(1) << (h & 63);
  return (tp.position_bits[static_cast<size_t>(pos) * tp.B + block] & bit) != 0;
}

// Exact O(P * number of affected support intervals) swap evaluation for equal pik.
double equal_swap_value(const EqualState& st,
                        const BalanceProblem& pr,
                        const EqualTemplate& tp,
                        int i,
                        int j,
                        std::vector<double>& d,
                        double& d2) {
  const int ui = st.sigma[i];
  const int uj = st.sigma[j];
  const double w = pr.inv_pi[0];
  d2 = 0.0;
  const size_t bi = static_cast<size_t>(ui) * pr.P;
  const size_t bj = static_cast<size_t>(uj) * pr.P;
  for (int q = 0; q < pr.P; ++q) {
    d[q] = w * (pr.X[bj + q] - pr.X[bi + q]) / pr.scale[q];
    d2 += d[q] * d[q];
  }

  double delta_f = 0.0;
  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0;
      uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const bool i_in = ((wi >> bit_index) & 1ULL) != 0ULL;
        const double sign = i_in ? 1.0 : -1.0;
        double dot = 0.0;
        const size_t base = static_cast<size_t>(h) * pr.P;
        for (int q = 0; q < pr.P; ++q) dot += st.residual[base + q] * d[q];
        delta_f += tp.probs[h] * (2.0 * sign * dot + d2);
      }
      x &= (x - 1ULL);
    }
  }
  return st.f + delta_f;
}

void apply_equal_swap(EqualState& st,
                      const BalanceProblem& pr,
                      const EqualTemplate& tp,
                      int i,
                      int j,
                      double new_f,
                      const std::vector<double>& d,
                      double d2) {

  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0;
      uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const bool i_in = ((wi >> bit_index) & 1ULL) != 0ULL;
        const double sign = i_in ? 1.0 : -1.0;
        const size_t base = static_cast<size_t>(h) * pr.P;
        double dot = 0.0;
        for (int q = 0; q < pr.P; ++q) dot += st.residual[base + q] * d[q];
        st.score[h] += 2.0 * sign * dot + d2;
        for (int q = 0; q < pr.P; ++q) st.residual[base + q] += sign * d[q];
      }
      x &= (x - 1ULL);
    }
  }
  std::swap(st.sigma[i], st.sigma[j]);
  st.f = new_f;
}

struct RunDiagnostics {
  long long adjacent_evals = 0;
  long long full_evals = 0;
  int adjacent_moves = 0;
  int global_moves = 0;
  int type1_hits = 0;
  bool type1_certified = false;
  bool type2_certified = false;
  bool timed_out = false;
};

inline bool is_circular_adjacent_pair(int i, int j, int N) {
  return (j == i + 1) || (i == 0 && j == N - 1);
}

template <typename CandidateFn, typename ApplyFn>
bool scan_adjacent_first(int N,
                         CandidateFn candidate,
                         ApplyFn apply,
                         double current_f,
                         double eps,
                         const std::chrono::steady_clock::time_point& t0,
                         double max_time,
                         RunDiagnostics& dg,
                         bool& interrupted) {
  interrupted = false;
  if (N <= 1) return false;
  long long local_count = 0;
  for (int i = 0; i <= N - 2; ++i) {
    const int j = i + 1;
    const double cf = candidate(i, j);
    ++dg.adjacent_evals; ++local_count;
    if (cf < current_f - eps) { apply(i, j, cf); return true; }
    if ((local_count & 255LL) == 0LL && std::isfinite(max_time)) {
      const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      if (e >= max_time) { interrupted = true; dg.timed_out = true; return false; }
    }
  }
  // Circular edge joining the last and first positions. For N=2 this is the
  // same transposition already tested above, so do not evaluate it twice.
  if (N > 2) {
    const double cf = candidate(N - 1, 0);
    ++dg.adjacent_evals;
    if (cf < current_f - eps) { apply(N - 1, 0, cf); return true; }
  }
  if (std::isfinite(max_time)) {
    const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (e >= max_time) { interrupted = true; dg.timed_out = true; }
  }
  return false;
}

template <typename CandidateFn, typename ApplyFn>
bool scan_full_first(int N,
                     CandidateFn candidate,
                     ApplyFn apply,
                     double current_f,
                     double eps,
                     bool skip_adjacent,
                     const std::chrono::steady_clock::time_point& t0,
                     double max_time,
                     RunDiagnostics& dg,
                     bool& interrupted) {
  interrupted = false;
  if (N <= 1) return false;
  long long local_count = 0;
  for (int i = 0; i <= N - 2; ++i) {
    for (int j = i + 1; j <= N - 1; ++j) {
      if (skip_adjacent && is_circular_adjacent_pair(i, j, N)) continue;
      const double cf = candidate(i, j);
      ++dg.full_evals; ++local_count;
      if (cf < current_f - eps) { apply(i, j, cf); return true; }
      if ((local_count & 255LL) == 0LL && std::isfinite(max_time)) {
        const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (e >= max_time) { interrupted = true; dg.timed_out = true; return false; }
      }
    }
  }
  if (std::isfinite(max_time)) {
    const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (e >= max_time) { interrupted = true; dg.timed_out = true; }
  }
  return false;
}

List run_balance_engine(std::vector<int> sigma,
                        const BalanceProblem& pr,
                        std::string search,
                        double eps,
                        double tol_pi,
                        int max_iter,
                        double max_time) {
  if (search != "adaptive" && search != "type1" && search != "type2") {
    stop("search must be one of 'adaptive', 'type1', or 'type2'.");
  }
  if (eps < 0.0) stop("improvement_tol must be non-negative.");
  if (max_iter < 0) stop("max_iter must be non-negative.");
  if (!(max_time > 0.0) && std::isfinite(max_time)) stop("max_time must be positive or Inf.");

  const bool equal_pi = all_equal_pi_std(pr.pi, tol_pi);
  RunDiagnostics dg;
  std::vector<double> f_path;
  std::vector<std::string> move_type;
  const auto t0 = std::chrono::steady_clock::now();

  int iter = 0;
  std::string stop_reason = "unknown";
  double f = std::numeric_limits<double>::quiet_NaN();
  double initial_f = std::numeric_limits<double>::quiet_NaN();

  // The equal-pik path carries the full state across accepted swaps. This avoids
  // rebuilding all systematic samples at every VNS iteration.
  if (equal_pi) {
    EqualTemplate tp = make_equal_template(pr, tol_pi);
    EqualState st = make_equal_state(sigma, pr, tp);
    // selected_positions is only needed to build the initial residual state.
    // Keep the compact bitset representation during the search.
    tp.selected_positions.clear();
    tp.selected_positions.shrink_to_fit();
    f = st.f;
    initial_f = f;
    f_path.push_back(f);

    std::vector<double> swap_delta(pr.P, 0.0);
    double swap_delta_sq = 0.0;
    auto candidate = [&](int i, int j) {
      return equal_swap_value(st, pr, tp, i, j, swap_delta, swap_delta_sq);
    };
    auto apply = [&](int i, int j, double nf) {
      apply_equal_swap(st, pr, tp, i, j, nf, swap_delta, swap_delta_sq);
    };

    while (true) {
      if (iter >= max_iter) { stop_reason = "max_iter"; break; }
      if (std::isfinite(max_time)) {
        const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (e >= max_time) { stop_reason = "max_time"; break; }
      }

      bool interrupted = false;
      if (search == "type2") {
        const bool moved = scan_full_first(pr.N, candidate, apply, st.f, eps, false, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason = "max_time"; break; }
        if (!moved) {
          dg.type2_certified = true;
          dg.type1_certified = true;
          stop_reason = "type2_local_minimum";
          break;
        }
        ++iter; ++dg.global_moves;
        dg.type1_certified = false;
        f_path.push_back(st.f); move_type.push_back("global");
        continue;
      }

      const bool moved_adj = scan_adjacent_first(pr.N, candidate, apply, st.f, eps, t0, max_time, dg, interrupted);
      if (interrupted) { stop_reason = "max_time"; break; }
      if (moved_adj) {
        ++iter; ++dg.adjacent_moves;
        dg.type1_certified = false;
        f_path.push_back(st.f); move_type.push_back("adjacent");
        continue;
      }

      dg.type1_certified = true;
      ++dg.type1_hits;
      if (search == "type1") {
        stop_reason = "type1_local_minimum";
        break;
      }

      const bool moved_full = scan_full_first(pr.N, candidate, apply, st.f, eps, true, t0, max_time, dg, interrupted);
      if (interrupted) { stop_reason = "max_time"; break; }
      if (!moved_full) {
        dg.type2_certified = true;
        stop_reason = "type2_local_minimum";
        break;
      }
      ++iter; ++dg.global_moves;
      dg.type1_certified = false;
      f_path.push_back(st.f); move_type.push_back("global");
    }
    sigma = st.sigma;
    f = st.f;
  } else {
    f = evaluate_balance_general(sigma, pr, tol_pi);
    initial_f = f;
    f_path.push_back(f);

    auto candidate = [&](int i, int j) {
      std::swap(sigma[i], sigma[j]);
      const double cf = evaluate_balance_general(sigma, pr, tol_pi);
      std::swap(sigma[i], sigma[j]);
      return cf;
    };
    auto apply = [&](int i, int j, double nf) { std::swap(sigma[i], sigma[j]); f = nf; };

    while (true) {
      if (iter >= max_iter) { stop_reason = "max_iter"; break; }
      if (std::isfinite(max_time)) {
        const double e = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (e >= max_time) { stop_reason = "max_time"; break; }
      }

      bool interrupted = false;
      if (search == "type2") {
        const bool moved = scan_full_first(pr.N, candidate, apply, f, eps, false, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason = "max_time"; break; }
        if (!moved) {
          dg.type2_certified = true;
          dg.type1_certified = true;
          stop_reason = "type2_local_minimum";
          break;
        }
        ++iter; ++dg.global_moves;
        dg.type1_certified = false;
        f_path.push_back(f); move_type.push_back("global");
        continue;
      }

      const bool moved_adj = scan_adjacent_first(pr.N, candidate, apply, f, eps, t0, max_time, dg, interrupted);
      if (interrupted) { stop_reason = "max_time"; break; }
      if (moved_adj) {
        ++iter; ++dg.adjacent_moves;
        dg.type1_certified = false;
        f_path.push_back(f); move_type.push_back("adjacent");
        continue;
      }

      dg.type1_certified = true;
      ++dg.type1_hits;
      if (search == "type1") {
        stop_reason = "type1_local_minimum";
        break;
      }

      const bool moved_full = scan_full_first(pr.N, candidate, apply, f, eps, true, t0, max_time, dg, interrupted);
      if (interrupted) { stop_reason = "max_time"; break; }
      if (!moved_full) {
        dg.type2_certified = true;
        stop_reason = "type2_local_minimum";
        break;
      }
      ++iter; ++dg.global_moves;
      dg.type1_certified = false;
      f_path.push_back(f); move_type.push_back("global");
    }
  }

  // Canonicalize the circular representation only after the search.  During
  // optimization all N positions are allowed to move; this is necessary if a
  // Type-II certificate is to mean all pairwise transpositions on the circle.
  const std::vector<int>::iterator anchor_it = std::find(sigma.begin(), sigma.end(), 0);
  if (anchor_it != sigma.end() && anchor_it != sigma.begin()) {
    std::rotate(sigma.begin(), anchor_it, sigma.end());
  }

  // One exact evaluation at the end guards against accumulated floating-point drift.
  const double f_exact = evaluate_balance_general(sigma, pr, tol_pi);
  f = f_exact;
  if (!f_path.empty()) f_path.back() = f_exact;

  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  IntegerVector sigma_out(pr.N);
  for (int i = 0; i < pr.N; ++i) sigma_out[i] = sigma[i] + 1;

  CharacterVector moves(move_type.size());
  for (size_t i = 0; i < move_type.size(); ++i) moves[i] = move_type[i];

  std::string local_type = "none";
  if (dg.type2_certified) local_type = "type2";
  else if (dg.type1_certified) local_type = "type1";

  return List::create(
    _["sigma"] = sigma_out,
    _["f"] = f,
    _["initial_f"] = initial_f,
    _["iterations"] = iter,
    _["elapsed_time"] = elapsed,
    _["f_path"] = f_path,
    _["move_type"] = moves,
    _["equal_pik"] = equal_pi,
    _["evaluation_method"] = equal_pi ? "incremental_equal_pik" : "general_unequal_pik",
    _["type1_certified"] = dg.type1_certified,
    _["type2_certified"] = dg.type2_certified,
    _["local_minimum_type"] = local_type,
    _["stop_reason"] = stop_reason,
    _["n_adjacent_evaluations"] = static_cast<double>(dg.adjacent_evals),
    _["n_full_evaluations"] = static_cast<double>(dg.full_evals),
    _["n_adjacent_moves"] = dg.adjacent_moves,
    _["n_global_moves"] = dg.global_moves,
    _["n_type1_certificates"] = dg.type1_hits
  );
}

// ============================================================
// Spatial spreading objectives
// ============================================================

struct SBProblem {
  int N;
  int K;
  int n;
  std::vector<double> pi;
  std::vector<double> coord;   // row-major N x K
  std::vector<double> dist2;   // optional dense squared-distance matrix
  bool has_dist2;
};

SBProblem make_sb_problem(const NumericVector& piR,
                          const NumericMatrix& coordR,
                          int dense_distance_limit = 2000) {
  const int N = piR.size();
  if (coordR.nrow() != N) stop("For objective='spread' and metric='sb', coordinates must have length(pik) rows.");
  const int K = coordR.ncol();
  if (K < 1) stop("coordinates must have at least one column.");

  SBProblem pr;
  pr.N = N;
  pr.K = K;
  pr.pi.assign(piR.begin(), piR.end());
  pr.coord.resize(static_cast<size_t>(N) * K);
  long double sump = 0.0L;
  for (int i = 0; i < N; ++i) {
    const double p = pr.pi[i];
    if (!std::isfinite(p) || p <= 0.0 || p > 1.0 + 1e-12) stop("All inclusion probabilities must be in (0,1].");
    sump += static_cast<long double>(p);
    for (int q = 0; q < K; ++q) {
      const double z = coordR(i, q);
      if (!std::isfinite(z)) stop("coordinates must contain only finite values.");
      pr.coord[static_cast<size_t>(i) * K + q] = z;
    }
  }
  if (std::fabs(static_cast<double>(sump - std::round(sump))) > 1e-8) {
    stop("sum(pik) must be an integer for fixed-size systematic sampling.");
  }
  pr.n = static_cast<int>(std::llround(static_cast<double>(sump)));
  if (pr.n <= 0 || pr.n >= N) stop("Spreading objectives require 0 < sum(pik) < N.");

  pr.has_dist2 = (N <= dense_distance_limit);
  if (pr.has_dist2) {
    pr.dist2.assign(static_cast<size_t>(N) * N, 0.0);
    for (int i = 0; i < N; ++i) {
      for (int j = i + 1; j < N; ++j) {
        double d = 0.0;
        const size_t bi = static_cast<size_t>(i) * K;
        const size_t bj = static_cast<size_t>(j) * K;
        for (int q = 0; q < K; ++q) {
          const double e = pr.coord[bi + q] - pr.coord[bj + q];
          d += e * e;
        }
        pr.dist2[static_cast<size_t>(i) * N + j] = d;
        pr.dist2[static_cast<size_t>(j) * N + i] = d;
      }
    }
  }
  return pr;
}

inline double sb_dist2(const SBProblem& pr, int i, int j) {
  if (pr.has_dist2) return pr.dist2[static_cast<size_t>(i) * pr.N + j];
  double d = 0.0;
  const size_t bi = static_cast<size_t>(i) * pr.K;
  const size_t bj = static_cast<size_t>(j) * pr.K;
  for (int q = 0; q < pr.K; ++q) {
    const double e = pr.coord[bi + q] - pr.coord[bj + q];
    d += e * e;
  }
  return d;
}

struct SBSupportState {
  std::vector<int> nearest;
  std::vector<int> second;
  std::vector<double> d1;
  std::vector<double> d2;
  std::vector<double> mass;
  double score;
};

SBSupportState build_sb_support_state(const std::vector<int>& selected_units,
                                      const SBProblem& pr) {
  if (selected_units.empty()) stop("A systematic support sample is empty.");
  SBSupportState st;
  st.nearest.assign(pr.N, -1);
  st.second.assign(pr.N, -1);
  st.d1.assign(pr.N, std::numeric_limits<double>::infinity());
  st.d2.assign(pr.N, std::numeric_limits<double>::infinity());
  st.mass.assign(pr.N, 0.0);

  for (int l = 0; l < pr.N; ++l) {
    int best = -1, second = -1;
    double dbest = std::numeric_limits<double>::infinity();
    double dsecond = std::numeric_limits<double>::infinity();
    for (int u : selected_units) {
      const double d = sb_dist2(pr, l, u);
      if (d < dbest) {
        dsecond = dbest; second = best;
        dbest = d; best = u;
      } else if (d < dsecond) {
        dsecond = d; second = u;
      }
    }
    st.nearest[l] = best;
    st.second[l] = second;
    st.d1[l] = dbest;
    st.d2[l] = dsecond;
    st.mass[best] += pr.pi[l];
  }

  double ss = 0.0;
  for (int u : selected_units) {
    const double e = st.mass[u] - 1.0;
    ss += e * e;
  }
  st.score = ss / static_cast<double>(pr.n);
  return st;
}

inline std::vector<int> units_from_positions(const std::vector<int>& sigma,
                                             const std::vector<int>& pos) {
  std::vector<int> units;
  units.reserve(pos.size());
  for (int p : pos) units.push_back(sigma[p]);
  return units;
}

double sb_sample_score(const std::vector<int>& selected_units,
                       const SBProblem& pr) {
  return build_sb_support_state(selected_units, pr).score;
}

// Evaluate one replacement a(out) -> b(in) for a fixed support sample.  The
// current nearest and second-nearest sampled units make this O(N) rather than
// O(N*n).  Ties have probability zero for ordinary continuous coordinates;
// the strict comparison gives deterministic behavior if ties occur.
double sb_replacement_score(const SBSupportState& st,
                            const SBProblem& pr,
                            int out_unit,
                            int in_unit,
                            std::vector<double>& delta_mass,
                            std::vector<unsigned char>& touched_flag,
                            std::vector<int>& touched) {
  touched.clear();
  auto add_delta = [&](int u, double z) {
    if (u < 0) return;
    if (!touched_flag[u]) {
      touched_flag[u] = 1;
      touched.push_back(u);
    }
    delta_mass[u] += z;
  };

  for (int l = 0; l < pr.N; ++l) {
    const int oldc = st.nearest[l];
    int basec = oldc;
    double based = st.d1[l];
    if (oldc == out_unit) {
      basec = st.second[l];
      based = st.d2[l];
    }
    const double din = sb_dist2(pr, l, in_unit);
    const int newc = (basec < 0 || din < based) ? in_unit : basec;
    if (newc != oldc) {
      add_delta(oldc, -pr.pi[l]);
      add_delta(newc, pr.pi[l]);
    }
  }

  double delta_ss = 0.0;
  for (int u : touched) {
    const double old_mass = st.mass[u];
    const double new_mass = old_mass + delta_mass[u];
    const bool old_selected = (u != in_unit);
    const bool new_selected = (u != out_unit);
    if (old_selected) {
      const double e = old_mass - 1.0;
      delta_ss -= e * e;
    }
    if (new_selected) {
      const double e = new_mass - 1.0;
      delta_ss += e * e;
    }
    delta_mass[u] = 0.0;
    touched_flag[u] = 0;
  }
  return st.score + delta_ss / static_cast<double>(pr.n);
}

struct SBEqualState {
  std::vector<int> sigma;
  std::vector<SBSupportState> support;
  double f;
};

SBEqualState make_sb_equal_state(const std::vector<int>& sigma,
                                 const SBProblem& pr,
                                 const EqualTemplate& tp) {
  SBEqualState st;
  st.sigma = sigma;
  st.support.resize(tp.D);
  st.f = 0.0;
  for (int h = 0; h < tp.D; ++h) {
    const std::vector<int> units = units_from_positions(sigma, tp.selected_positions[h]);
    st.support[h] = build_sb_support_state(units, pr);
    st.f += tp.probs[h] * st.support[h].score;
  }
  return st;
}

double sb_equal_swap_value(const SBEqualState& st,
                           const SBProblem& pr,
                           const EqualTemplate& tp,
                           int i,
                           int j,
                           std::vector<double>& delta_mass,
                           std::vector<unsigned char>& touched_flag,
                           std::vector<int>& touched) {
  const int ui = st.sigma[i];
  const int uj = st.sigma[j];
  double nf = st.f;
  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0; uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const bool i_in = ((wi >> bit_index) & 1ULL) != 0ULL;
        const int out_unit = i_in ? ui : uj;
        const int in_unit  = i_in ? uj : ui;
        const double sc = sb_replacement_score(st.support[h], pr, out_unit, in_unit,
                                               delta_mass, touched_flag, touched);
        nf += tp.probs[h] * (sc - st.support[h].score);
      }
      x &= (x - 1ULL);
    }
  }
  return nf;
}

void apply_sb_equal_swap(SBEqualState& st,
                         const SBProblem& pr,
                         const EqualTemplate& tp,
                         int i,
                         int j,
                         double /*new_f*/) {
  std::swap(st.sigma[i], st.sigma[j]);
  double fnew = st.f;
  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0; uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const double old_sc = st.support[h].score;
        const std::vector<int> units = units_from_positions(st.sigma, tp.selected_positions[h]);
        st.support[h] = build_sb_support_state(units, pr);
        fnew += tp.probs[h] * (st.support[h].score - old_sc);
      }
      x &= (x - 1ULL);
    }
  }
  st.f = fnew;
}

struct WEdge { int idx; double w; };

struct SparseW {
  int N;
  std::vector<std::vector<WEdge>> row;
  std::vector<std::vector<WEdge>> col;
  std::vector<double> rowsum;
  std::vector<double> colsum;
  double total;
};

SparseW make_sparse_w(int N,
                      const IntegerVector& wiR,
                      const IntegerVector& wjR,
                      const NumericVector& wxR) {
  if (wiR.size() != wjR.size() || wiR.size() != wxR.size()) stop("Invalid sparse W triplets.");
  SparseW W;
  W.N = N;
  W.row.resize(N);
  W.col.resize(N);
  W.rowsum.assign(N, 0.0);
  W.colsum.assign(N, 0.0);
  W.total = 0.0;
  for (R_xlen_t k = 0; k < wxR.size(); ++k) {
    const int i = wiR[k] - 1;
    const int j = wjR[k] - 1;
    const double w = wxR[k];
    if (i < 0 || i >= N || j < 0 || j >= N || !std::isfinite(w)) stop("Invalid sparse W entry.");
    if (w == 0.0) continue;
    W.row[i].push_back(WEdge{j, w});
    W.col[j].push_back(WEdge{i, w});
    W.rowsum[i] += w;
    W.colsum[j] += w;
    W.total += w;
  }
  if (!(W.total > 0.0)) stop("W must have positive total weight.");
  for (int i = 0; i < N; ++i) {
    if (!(W.rowsum[i] > 0.0)) stop("Every row of W must have positive weight for the IB criterion.");
    std::sort(W.row[i].begin(), W.row[i].end(), [](const WEdge& a, const WEdge& b){ return a.idx < b.idx; });
    std::sort(W.col[i].begin(), W.col[i].end(), [](const WEdge& a, const WEdge& b){ return a.idx < b.idx; });
  }
  return W;
}

inline double get_w(const SparseW& W, int i, int j) {
  const std::vector<WEdge>& r = W.row[i];
  for (const WEdge& e : r) {
    if (e.idx == j) return e.w;
    if (e.idx > j) break;
  }
  return 0.0;
}

inline double ib_from_scalars(double A0, double R, double C, double Q, const SparseW& W) {
  const double T = W.total;
  const double num = A0 - (R * C) / T;
  const double d1 = R - (R * R) / T;
  const double d2 = Q - (C * C) / T;
  const double den2 = d1 * d2;
  if (!(den2 > 1e-18) || !std::isfinite(den2)) return std::numeric_limits<double>::infinity();
  return num / std::sqrt(den2);
}

struct IBSupportState {
  std::vector<double> ws;   // W s
  std::vector<double> wts;  // W' s
  double A0;
  double R;
  double C;
  double Q;
  double score;
};

IBSupportState build_ib_support_state(const std::vector<int>& selected_units,
                                      const SparseW& W) {
  IBSupportState st;
  st.ws.assign(W.N, 0.0);
  st.wts.assign(W.N, 0.0);
  st.R = 0.0; st.C = 0.0;
  std::vector<unsigned char> selected(W.N, 0);
  for (int u : selected_units) {
    selected[u] = 1;
    st.R += W.rowsum[u];
    st.C += W.colsum[u];
    for (const WEdge& e : W.col[u]) st.ws[e.idx] += e.w;
    for (const WEdge& e : W.row[u]) st.wts[e.idx] += e.w;
  }
  st.A0 = 0.0;
  for (int u : selected_units) st.A0 += st.ws[u];
  st.Q = 0.0;
  for (int i = 0; i < W.N; ++i) st.Q += st.ws[i] * st.ws[i] / W.rowsum[i];
  st.score = ib_from_scalars(st.A0, st.R, st.C, st.Q, W);
  return st;
}

double ib_sample_score(const std::vector<int>& selected_units,
                       const SparseW& W) {
  return build_ib_support_state(selected_units, W).score;
}

struct IBCandidate {
  double A0;
  double R;
  double C;
  double Q;
  double score;
};

IBCandidate ib_replacement_candidate(const IBSupportState& st,
                                     const SparseW& W,
                                     int out_unit,
                                     int in_unit) {
  IBCandidate z;
  z.R = st.R - W.rowsum[out_unit] + W.rowsum[in_unit];
  z.C = st.C - W.colsum[out_unit] + W.colsum[in_unit];
  const double cross = get_w(W, in_unit, in_unit) - get_w(W, in_unit, out_unit)
                     - get_w(W, out_unit, in_unit) + get_w(W, out_unit, out_unit);
  z.A0 = st.A0 + (st.ws[in_unit] - st.ws[out_unit])
                 + (st.wts[in_unit] - st.wts[out_unit]) + cross;

  z.Q = st.Q;
  const std::vector<WEdge>& ca = W.col[out_unit];
  const std::vector<WEdge>& cb = W.col[in_unit];
  size_t ia = 0, ib = 0;
  while (ia < ca.size() || ib < cb.size()) {
    int row;
    double wa = 0.0, wb = 0.0;
    if (ib >= cb.size() || (ia < ca.size() && ca[ia].idx < cb[ib].idx)) {
      row = ca[ia].idx; wa = ca[ia].w; ++ia;
    } else if (ia >= ca.size() || cb[ib].idx < ca[ia].idx) {
      row = cb[ib].idx; wb = cb[ib].w; ++ib;
    } else {
      row = ca[ia].idx; wa = ca[ia].w; wb = cb[ib].w; ++ia; ++ib;
    }
    const double d = wb - wa;
    z.Q += (2.0 * st.ws[row] * d + d * d) / W.rowsum[row];
  }
  z.score = ib_from_scalars(z.A0, z.R, z.C, z.Q, W);
  return z;
}

struct IBEqualState {
  std::vector<int> sigma;
  std::vector<IBSupportState> support;
  double f;
};

IBEqualState make_ib_equal_state(const std::vector<int>& sigma,
                                 const SparseW& W,
                                 const EqualTemplate& tp) {
  IBEqualState st;
  st.sigma = sigma;
  st.support.resize(tp.D);
  st.f = 0.0;
  for (int h = 0; h < tp.D; ++h) {
    const std::vector<int> units = units_from_positions(sigma, tp.selected_positions[h]);
    st.support[h] = build_ib_support_state(units, W);
    st.f += tp.probs[h] * st.support[h].score;
  }
  return st;
}

double ib_equal_swap_value(const IBEqualState& st,
                           const SparseW& W,
                           const EqualTemplate& tp,
                           int i,
                           int j,
                           std::vector<IBCandidate>& scratch) {
  const int ui = st.sigma[i];
  const int uj = st.sigma[j];
  double nf = st.f;
  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0; uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const bool i_in = ((wi >> bit_index) & 1ULL) != 0ULL;
        const int out_unit = i_in ? ui : uj;
        const int in_unit  = i_in ? uj : ui;
        scratch[h] = ib_replacement_candidate(st.support[h], W, out_unit, in_unit);
        nf += tp.probs[h] * (scratch[h].score - st.support[h].score);
      }
      x &= (x - 1ULL);
    }
  }
  return nf;
}

void apply_ib_equal_swap(IBEqualState& st,
                         const SparseW& W,
                         const EqualTemplate& tp,
                         int i,
                         int j,
                         double new_f,
                         const std::vector<IBCandidate>& scratch) {
  const int ui = st.sigma[i];
  const int uj = st.sigma[j];
  for (int block = 0; block < tp.B; ++block) {
    const uint64_t wi = tp.position_bits[static_cast<size_t>(i) * tp.B + block];
    const uint64_t wj = tp.position_bits[static_cast<size_t>(j) * tp.B + block];
    uint64_t x = wi ^ wj;
    while (x) {
#if defined(__GNUG__)
      const int bit_index = __builtin_ctzll(x);
#else
      int bit_index = 0; uint64_t y = x;
      while ((y & 1ULL) == 0ULL) { y >>= 1; ++bit_index; }
#endif
      const int h = (block << 6) + bit_index;
      if (h < tp.D) {
        const bool i_in = ((wi >> bit_index) & 1ULL) != 0ULL;
        const int out_unit = i_in ? ui : uj;
        const int in_unit  = i_in ? uj : ui;
        IBSupportState& ss = st.support[h];
        for (const WEdge& e : W.col[out_unit]) ss.ws[e.idx] -= e.w;
        for (const WEdge& e : W.col[in_unit])  ss.ws[e.idx] += e.w;
        for (const WEdge& e : W.row[out_unit]) ss.wts[e.idx] -= e.w;
        for (const WEdge& e : W.row[in_unit])  ss.wts[e.idx] += e.w;
        ss.A0 = scratch[h].A0;
        ss.R = scratch[h].R;
        ss.C = scratch[h].C;
        ss.Q = scratch[h].Q;
        ss.score = scratch[h].score;
      }
      x &= (x - 1ULL);
    }
  }
  std::swap(st.sigma[i], st.sigma[j]);
  st.f = new_f;
}

std::vector<double> systematic_cumV(const std::vector<int>& sigma,
                                    const std::vector<double>& pi,
                                    bool equal_pi) {
  const int N = static_cast<int>(pi.size());
  std::vector<double> cumV(N);
  if (equal_pi) {
    const double p0 = pi[0];
    for (int k = 0; k < N; ++k) cumV[k] = static_cast<double>(k + 1) * p0;
  } else {
    long double acc = 0.0L;
    for (int k = 0; k < N; ++k) { acc += static_cast<long double>(pi[sigma[k]]); cumV[k] = static_cast<double>(acc); }
  }
  return cumV;
}

double evaluate_spread_general(const std::vector<int>& sigma,
                               const std::vector<double>& pi,
                               const std::string& metric,
                               const SBProblem* sb,
                               const SparseW* W,
                               double tol_pi) {
  const bool equal_pi = all_equal_pi_std(pi, tol_pi);
  const std::vector<double> cumV = systematic_cumV(sigma, pi, equal_pi);
  std::vector<double> cuts;
  if (equal_pi) {
    const int N = static_cast<int>(pi.size());
    const int n = static_cast<int>(std::llround(pi[0] * static_cast<double>(N)));
    const int D = N / std::max(gcd_int(N, n), 1);
    cuts.resize(D + 1);
    for (int h = 0; h <= D; ++h) cuts[h] = static_cast<double>(h) / static_cast<double>(D);
  } else {
    cuts = build_cuts(cumV, tol_pi);
  }

  double f = 0.0;
  for (size_t h = 0; h + 1 < cuts.size(); ++h) {
    const double prob = cuts[h + 1] - cuts[h];
    const double u = 0.5 * (cuts[h + 1] + cuts[h]);
    const std::vector<int> pos = selected_positions_from_u(cumV, u);
    const std::vector<int> units = units_from_positions(sigma, pos);
    const double sc = (metric == "sb") ? sb_sample_score(units, *sb) : ib_sample_score(units, *W);
    f += prob * sc;
  }
  return f;
}

List run_spread_engine(std::vector<int> sigma,
                       const std::vector<double>& pi,
                       const std::string& metric,
                       const SBProblem* sb,
                       const SparseW* W,
                       std::string search,
                       double eps,
                       double tol_pi,
                       int max_iter,
                       double max_time,
                       double max_state_cells) {
  if (search != "adaptive" && search != "type1" && search != "type2") stop("search must be one of 'adaptive', 'type1', or 'type2'.");
  if (metric != "sb" && metric != "ib") stop("spread metric must be 'sb' or 'ib'.");
  const int N = static_cast<int>(pi.size());
  const bool equal_pi = all_equal_pi_std(pi, tol_pi);
  bool use_equal_incremental = equal_pi;
  if (equal_pi) {
    const int n0 = static_cast<int>(std::llround(pi[0] * static_cast<double>(N)));
    const int D0 = N / std::max(gcd_int(N, n0), 1);
    const long double cells = static_cast<long double>(N) * static_cast<long double>(D0);
    if (std::isfinite(max_state_cells) && cells > max_state_cells) use_equal_incremental = false;
  }
  RunDiagnostics dg;
  std::vector<double> f_path;
  std::vector<std::string> move_type;
  const auto t0 = std::chrono::steady_clock::now();
  int iter = 0;
  std::string stop_reason = "unknown";
  double f = std::numeric_limits<double>::quiet_NaN();
  double initial_f = std::numeric_limits<double>::quiet_NaN();

  if (use_equal_incremental) {
    // Equal-pik systematic support is invariant to sigma. This lets both
    // spreading criteria update only support samples containing exactly one
    // of the swapped positions.
    BalanceProblem dummy;
    dummy.N = N; dummy.P = 1; dummy.pi = pi; dummy.inv_pi.resize(N);
    for (int i = 0; i < N; ++i) dummy.inv_pi[i] = 1.0 / pi[i];
    EqualTemplate tp = make_equal_template(dummy, tol_pi);

    if (metric == "sb") {
      SBEqualState st = make_sb_equal_state(sigma, *sb, tp);
      f = st.f; initial_f = f; f_path.push_back(f);
      std::vector<double> delta_mass(N, 0.0);
      std::vector<unsigned char> touched_flag(N, 0);
      std::vector<int> touched; touched.reserve(N);
      auto candidate = [&](int i, int j) {
        return sb_equal_swap_value(st, *sb, tp, i, j, delta_mass, touched_flag, touched);
      };
      auto apply = [&](int i, int j, double nf) { apply_sb_equal_swap(st, *sb, tp, i, j, nf); };

      while (true) {
        if (iter >= max_iter) { stop_reason = "max_iter"; break; }
        if (std::isfinite(max_time) && std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count() >= max_time) { stop_reason = "max_time"; break; }
        bool interrupted = false;
        if (search == "type2") {
          const bool moved = scan_full_first(N, candidate, apply, st.f, eps, false, t0, max_time, dg, interrupted);
          if (interrupted) { stop_reason = "max_time"; break; }
          if (!moved) { dg.type2_certified = true; dg.type1_certified = true; stop_reason = "type2_local_minimum"; break; }
          ++iter; ++dg.global_moves; dg.type1_certified = false; f_path.push_back(st.f); move_type.push_back("global"); continue;
        }
        const bool moved_adj = scan_adjacent_first(N, candidate, apply, st.f, eps, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason = "max_time"; break; }
        if (moved_adj) { ++iter; ++dg.adjacent_moves; dg.type1_certified=false; f_path.push_back(st.f); move_type.push_back("adjacent"); continue; }
        dg.type1_certified = true; ++dg.type1_hits;
        if (search == "type1") { stop_reason = "type1_local_minimum"; break; }
        const bool moved_full = scan_full_first(N, candidate, apply, st.f, eps, true, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason = "max_time"; break; }
        if (!moved_full) { dg.type2_certified=true; stop_reason="type2_local_minimum"; break; }
        ++iter; ++dg.global_moves; dg.type1_certified=false; f_path.push_back(st.f); move_type.push_back("global");
      }
      sigma = st.sigma; f = st.f;
    } else {
      IBEqualState st = make_ib_equal_state(sigma, *W, tp);
      f = st.f; initial_f = f; f_path.push_back(f);
      std::vector<IBCandidate> scratch(tp.D);
      auto candidate = [&](int i, int j) { return ib_equal_swap_value(st, *W, tp, i, j, scratch); };
      auto apply = [&](int i, int j, double nf) { apply_ib_equal_swap(st, *W, tp, i, j, nf, scratch); };

      while (true) {
        if (iter >= max_iter) { stop_reason = "max_iter"; break; }
        if (std::isfinite(max_time) && std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count() >= max_time) { stop_reason = "max_time"; break; }
        bool interrupted = false;
        if (search == "type2") {
          const bool moved = scan_full_first(N, candidate, apply, st.f, eps, false, t0, max_time, dg, interrupted);
          if (interrupted) { stop_reason = "max_time"; break; }
          if (!moved) { dg.type2_certified=true; dg.type1_certified=true; stop_reason="type2_local_minimum"; break; }
          ++iter; ++dg.global_moves; dg.type1_certified=false; f_path.push_back(st.f); move_type.push_back("global"); continue;
        }
        const bool moved_adj = scan_adjacent_first(N, candidate, apply, st.f, eps, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason = "max_time"; break; }
        if (moved_adj) { ++iter; ++dg.adjacent_moves; dg.type1_certified=false; f_path.push_back(st.f); move_type.push_back("adjacent"); continue; }
        dg.type1_certified=true; ++dg.type1_hits;
        if (search == "type1") { stop_reason="type1_local_minimum"; break; }
        const bool moved_full = scan_full_first(N, candidate, apply, st.f, eps, true, t0, max_time, dg, interrupted);
        if (interrupted) { stop_reason="max_time"; break; }
        if (!moved_full) { dg.type2_certified=true; stop_reason="type2_local_minimum"; break; }
        ++iter; ++dg.global_moves; dg.type1_certified=false; f_path.push_back(st.f); move_type.push_back("global");
      }
      sigma = st.sigma; f = st.f;
    }
  } else {
    f = evaluate_spread_general(sigma, pi, metric, sb, W, tol_pi);
    initial_f = f; f_path.push_back(f);
    auto candidate = [&](int i, int j) {
      std::swap(sigma[i], sigma[j]);
      const double cf = evaluate_spread_general(sigma, pi, metric, sb, W, tol_pi);
      std::swap(sigma[i], sigma[j]);
      return cf;
    };
    auto apply = [&](int i, int j, double nf) { std::swap(sigma[i], sigma[j]); f = nf; };
    while (true) {
      if (iter >= max_iter) { stop_reason="max_iter"; break; }
      if (std::isfinite(max_time) && std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count() >= max_time) { stop_reason="max_time"; break; }
      bool interrupted=false;
      if (search=="type2") {
        const bool moved=scan_full_first(N,candidate,apply,f,eps,false,t0,max_time,dg,interrupted);
        if (interrupted) { stop_reason="max_time"; break; }
        if (!moved) { dg.type2_certified=true; dg.type1_certified=true; stop_reason="type2_local_minimum"; break; }
        ++iter; ++dg.global_moves; dg.type1_certified=false; f_path.push_back(f); move_type.push_back("global"); continue;
      }
      const bool moved_adj=scan_adjacent_first(N,candidate,apply,f,eps,t0,max_time,dg,interrupted);
      if (interrupted) { stop_reason="max_time"; break; }
      if (moved_adj) { ++iter; ++dg.adjacent_moves; dg.type1_certified=false; f_path.push_back(f); move_type.push_back("adjacent"); continue; }
      dg.type1_certified=true; ++dg.type1_hits;
      if (search=="type1") { stop_reason="type1_local_minimum"; break; }
      const bool moved_full=scan_full_first(N,candidate,apply,f,eps,true,t0,max_time,dg,interrupted);
      if (interrupted) { stop_reason="max_time"; break; }
      if (!moved_full) { dg.type2_certified=true; stop_reason="type2_local_minimum"; break; }
      ++iter; ++dg.global_moves; dg.type1_certified=false; f_path.push_back(f); move_type.push_back("global");
    }
  }

  const std::vector<int>::iterator anchor_it = std::find(sigma.begin(), sigma.end(), 0);
  if (anchor_it != sigma.end() && anchor_it != sigma.begin()) std::rotate(sigma.begin(), anchor_it, sigma.end());
  const double f_exact = evaluate_spread_general(sigma, pi, metric, sb, W, tol_pi);
  f = f_exact;
  if (!f_path.empty()) f_path.back() = f_exact;

  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  IntegerVector sigma_out(N);
  for (int i=0;i<N;++i) sigma_out[i]=sigma[i]+1;
  CharacterVector moves(move_type.size());
  for (size_t i=0;i<move_type.size();++i) moves[i]=move_type[i];
  std::string local_type="none";
  if (dg.type2_certified) local_type="type2"; else if (dg.type1_certified) local_type="type1";
  std::string method;
  if (use_equal_incremental) method = (metric=="sb") ? "incremental_spread_sb_equal_pik" : "incremental_spread_ib_equal_pik";
  else if (equal_pi) method = "general_spread_equal_pik_memory_fallback";
  else method = "general_spread_unequal_pik";

  return List::create(
    _["sigma"]=sigma_out, _["f"]=f, _["initial_f"]=initial_f,
    _["iterations"]=iter, _["elapsed_time"]=elapsed, _["f_path"]=f_path,
    _["move_type"]=moves, _["equal_pik"]=equal_pi, _["evaluation_method"]=method,
    _["type1_certified"]=dg.type1_certified, _["type2_certified"]=dg.type2_certified,
    _["local_minimum_type"]=local_type, _["stop_reason"]=stop_reason,
    _["n_adjacent_evaluations"]=static_cast<double>(dg.adjacent_evals),
    _["n_full_evaluations"]=static_cast<double>(dg.full_evals),
    _["n_adjacent_moves"]=dg.adjacent_moves, _["n_global_moves"]=dg.global_moves,
    _["n_type1_certificates"]=dg.type1_hits
  );
}


} // namespace

// [[Rcpp::export]]
bool vns_all_equal_pik_cpp(NumericVector pik, double tol = 1e-12) {
  std::vector<double> p(pik.begin(), pik.end());
  return all_equal_pi_std(p, tol);
}

// [[Rcpp::export]]
double vns_balance_objective_cpp(IntegerVector sigma,
                                 NumericVector pik,
                                 NumericMatrix X,
                                 Nullable<NumericVector> balance_scale = R_NilValue,
                                 double tol_pik = 1e-12) {
  BalanceProblem pr = make_problem(pik, X, balance_scale);
  std::vector<int> s = validate_sigma(sigma, pr.N);
  return evaluate_balance_general(s, pr, tol_pik);
}

// [[Rcpp::export]]
List vns_run_balance_cpp(IntegerVector sigma,
                         NumericVector pik,
                         NumericMatrix X,
                         std::string search = "adaptive",
                         Nullable<NumericVector> balance_scale = R_NilValue,
                         double improvement_tol = 1e-10,
                         double tol_pik = 1e-12,
                         int max_iter = 100000,
                         double max_time = 1e300) {
  BalanceProblem pr = make_problem(pik, X, balance_scale);
  std::vector<int> s = validate_sigma(sigma, pr.N);
  return run_balance_engine(s, pr, search, improvement_tol, tol_pik, max_iter, max_time);
}


// [[Rcpp::export]]
double vns_spread_objective_cpp(IntegerVector sigma,
                                NumericVector pik,
                                NumericMatrix coordinates,
                                std::string metric,
                                IntegerVector wi,
                                IntegerVector wj,
                                NumericVector wx,
                                double tol_pik = 1e-12,
                                int dense_distance_limit = 2000) {
  std::vector<double> p(pik.begin(), pik.end());
  std::vector<int> s = validate_sigma(sigma, pik.size());
  if (metric == "sb") {
    SBProblem sb = make_sb_problem(pik, coordinates, dense_distance_limit);
    return evaluate_spread_general(s, p, metric, &sb, nullptr, tol_pik);
  }
  if (metric == "ib") {
    SparseW W = make_sparse_w(pik.size(), wi, wj, wx);
    return evaluate_spread_general(s, p, metric, nullptr, &W, tol_pik);
  }
  stop("spread metric must be 'sb' or 'ib'.");
  return NA_REAL;
}

// [[Rcpp::export]]
List vns_run_spread_cpp(IntegerVector sigma,
                        NumericVector pik,
                        NumericMatrix coordinates,
                        std::string metric,
                        IntegerVector wi,
                        IntegerVector wj,
                        NumericVector wx,
                        std::string search = "adaptive",
                        double improvement_tol = 1e-10,
                        double tol_pik = 1e-12,
                        int max_iter = 100000,
                        double max_time = 1e300,
                        int dense_distance_limit = 2000,
                        double max_state_cells = 5000000.0) {
  std::vector<double> p(pik.begin(), pik.end());
  std::vector<int> s = validate_sigma(sigma, pik.size());
  if (metric == "sb") {
    SBProblem sb = make_sb_problem(pik, coordinates, dense_distance_limit);
    return run_spread_engine(s, p, metric, &sb, nullptr, search, improvement_tol,
                             tol_pik, max_iter, max_time, max_state_cells);
  }
  if (metric == "ib") {
    SparseW W = make_sparse_w(pik.size(), wi, wj, wx);
    return run_spread_engine(s, p, metric, nullptr, &W, search, improvement_tol,
                             tol_pik, max_iter, max_time, max_state_cells);
  }
  stop("spread metric must be 'sb' or 'ib'.");
  return List::create();
}


// Exact support probabilities of a systematic design. Sample keys use 1-based unit labels.
// [[Rcpp::export]]
List systematic_support_cpp(IntegerVector sigma,
                            NumericVector pik,
                            double tol_pik = 1e-12) {
  const int N = pik.size();
  if (sigma.size() != N) stop("sigma and pik must have the same length.");
  std::vector<int> s = validate_sigma(sigma, N);
  std::vector<double> p(pik.begin(), pik.end());
  long double sump = 0.0L;
  for (double z : p) {
    if (!std::isfinite(z) || z <= 0.0 || z > 1.0 + 1e-12) stop("pik must lie in (0,1].");
    sump += static_cast<long double>(z);
  }
  if (std::fabs(static_cast<double>(sump - std::round(sump))) > 1e-8) stop("sum(pik) must be an integer.");

  const bool equal_pik = all_equal_pi_std(p, tol_pik);
  std::vector<double> cumV(N);
  std::vector<double> cuts;
  if (equal_pik) {
    const int n = static_cast<int>(std::llround(p[0] * static_cast<double>(N)));
    const int D0 = N / std::max(gcd_int(N, n), 1);
    const double p0 = static_cast<double>(n) / static_cast<double>(N);
    for (int k = 0; k < N; ++k) cumV[k] = static_cast<double>(k + 1) * p0;
    cuts.resize(D0 + 1);
    for (int h = 0; h <= D0; ++h) cuts[h] = static_cast<double>(h) / static_cast<double>(D0);
  } else {
    long double acc = 0.0L;
    for (int k = 0; k < N; ++k) { acc += static_cast<long double>(p[s[k]]); cumV[k] = static_cast<double>(acc); }
    cuts = build_cuts(cumV, tol_pik);
  }
  const int D = static_cast<int>(cuts.size()) - 1;
  List samples(D);
  NumericVector probs(D);

  for (int h = 0; h < D; ++h) {
    probs[h] = cuts[h + 1] - cuts[h];
    const double u = 0.5 * (cuts[h] + cuts[h + 1]);
    const std::vector<int> pos = selected_positions_from_u(cumV, u);
    IntegerVector units(pos.size());
    for (size_t q = 0; q < pos.size(); ++q) units[q] = s[pos[q]] + 1;
    std::sort(units.begin(), units.end());
    samples[h] = units;
  }
  return List::create(_["samples"] = samples, _["probabilities"] = probs);
}

// Exact second-order inclusion probabilities for a selected set of pairs.
// [[Rcpp::export]]
NumericVector systematic_pair_probs_cpp(IntegerVector sigma,
                                        NumericVector pik,
                                        IntegerMatrix pairs,
                                        double tol_pik = 1e-12) {
  const int N = pik.size();
  if (pairs.ncol() != 2) stop("pairs must be a two-column integer matrix.");
  std::vector<int> s = validate_sigma(sigma, N);
  std::vector<double> p(pik.begin(), pik.end());
  long double sump = 0.0L;
  for (double z : p) {
    if (!std::isfinite(z) || z <= 0.0 || z > 1.0 + 1e-12) stop("pik must lie in (0,1].");
    sump += static_cast<long double>(z);
  }
  if (std::fabs(static_cast<double>(sump - std::round(sump))) > 1e-8) stop("sum(pik) must be an integer.");

  const int K = pairs.nrow();
  std::vector<int> a(K), b(K);
  for (int q = 0; q < K; ++q) {
    a[q] = pairs(q,0) - 1; b[q] = pairs(q,1) - 1;
    if (a[q] < 0 || b[q] < 0 || a[q] >= N || b[q] >= N || a[q] == b[q]) stop("Invalid pair index.");
  }

  const bool equal_pik = all_equal_pi_std(p, tol_pik);
  std::vector<double> cumV(N);
  if (equal_pik) {
    const int n = static_cast<int>(std::llround(p[0] * static_cast<double>(N)));
    const double p0 = static_cast<double>(n) / static_cast<double>(N);
    for (int k = 0; k < N; ++k) cumV[k] = static_cast<double>(k + 1) * p0;
  } else {
    long double acc = 0.0L;
    for (int k = 0; k < N; ++k) {
      acc += static_cast<long double>(p[s[k]]);
      cumV[k] = static_cast<double>(acc);
    }
  }
  NumericVector out(K);
  std::vector<unsigned char> selected(N, 0);

  // Use the same breakpoint sweep as the general objective evaluator.  This
  // avoids rebuilding an O(N) membership vector independently on every
  // support interval, which is particularly important for unequal pik.
  const std::vector<BreakEvent> events = make_break_events(cumV, tol_pik);
  const double first_right = events.empty() ? 1.0 : events.front().r;
  const std::vector<int> initial_pos = selected_positions_from_u(cumV, 0.5 * first_right);
  for (int pos : initial_pos) selected[s[pos]] = 1;

  auto add_interval = [&](double prob) {
    if (prob <= 0.0) return;
    for (int q = 0; q < K; ++q) {
      if (selected[a[q]] && selected[b[q]]) out[q] += prob;
    }
  };

  double left = 0.0;
  size_t idx = 0;
  while (idx < events.size()) {
    const double r = events[idx].r;
    add_interval(r - left);
    size_t jdx = idx;
    while (jdx < events.size() && std::fabs(events[jdx].r - r) <= tol_pik) ++jdx;

    // As in the balance evaluator, collapse each consecutive chain of events
    // to its net replacement.  Internal pi=1 units remain selected throughout.
    size_t qevent = idx;
    while (qevent < jdx) {
      int a = events[qevent].k;
      int b = a;
      ++qevent;
      while (qevent < jdx && events[qevent].k == b + 1) {
        b = events[qevent].k;
        ++qevent;
      }
      selected[s[a]] = 0;
      selected[s[b + 1]] = 1;
    }
    left = r;
    idx = jdx;
  }
  add_interval(1.0 - left);
  return out;
}

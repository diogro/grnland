// Hill-type gene regulatory networks: development to steady state and
// analytical derivatives.
//
// Dynamics, for genes i = 1..n:
//   dx_i/dt = alpha_i * prod_j H_ij(x_j) - gamma_i * x_i
// with
//   H_ij(x) = a(x)      if regulation(i, j) == +1  (j activates i)
//   H_ij(x) = 1 - a(x)  if regulation(i, j) == -1  (j represses i)
//   H_ij(x) = 1         if regulation(i, j) ==  0
//   a(x)    = x^h / (theta^h + x^h), with theta = theta(i, j), h = hill(i, j).
//
// Loci are evolvable parameters, each identified by (code, i, j):
//   code 1 = alpha[i], 2 = gamma[i], 3 = theta[i, j], 4 = hill[i, j].
// Indices passed from R are 0-based.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

struct Net {
  int n;
  arma::vec alpha, gamma;
  arma::mat theta, hill;
  arma::imat reg;
};

Net make_net(const List& base) {
  Net m;
  m.alpha = as<arma::vec>(base["alpha"]);
  m.gamma = as<arma::vec>(base["gamma"]);
  m.theta = as<arma::mat>(base["theta"]);
  m.hill = as<arma::mat>(base["hill"]);
  m.reg = as<arma::imat>(base["regulation"]);
  m.n = m.alpha.n_elem;
  return m;
}

struct Loci {
  IntegerVector code, i, j;
  int size() const { return code.size(); }
};

void apply_genotype(Net& m, const Loci& loci, const double* g, int stride) {
  for (int l = 0; l < loci.size(); ++l) {
    double v = g[l * stride];
    int i = loci.i[l], j = loci.j[l];
    switch (loci.code[l]) {
    case 1: m.alpha(i) = v; break;
    case 2: m.gamma(i) = v; break;
    case 3: m.theta(i, j) = v; break;
    case 4: m.hill(i, j) = v; break;
    default: Rcpp::stop("unknown locus code");
    }
  }
}

// Activation function a(x) and its partial derivatives.
inline double act(double x, double th, double h) {
  if (x <= 0) return 0.0;
  double xn = std::pow(x, h), tn = std::pow(th, h);
  return xn / (tn + xn);
}
inline double dact_dx(double x, double th, double h) {
  if (x < 0) x = 0;
  double tn = std::pow(th, h), xn = std::pow(x, h);
  double d = tn + xn;
  return h * tn * std::pow(x, h - 1.0) / (d * d);
}
inline double dact_dth(double x, double th, double h) {
  if (x <= 0) return 0.0;
  double tn = std::pow(th, h), xn = std::pow(x, h);
  double d = tn + xn;
  return -h * std::pow(th, h - 1.0) * xn / (d * d);
}
inline double dact_dh(double x, double th, double h) {
  if (x <= 0) return 0.0;
  double tn = std::pow(th, h), xn = std::pow(x, h);
  double d = tn + xn;
  return xn * tn * (std::log(x) - std::log(th)) / (d * d);
}

inline double H(const Net& m, int i, int j, double x) {
  int r = m.reg(i, j);
  if (r == 0) return 1.0;
  double a = act(x, m.theta(i, j), m.hill(i, j));
  return r > 0 ? a : 1.0 - a;
}

// Product of H_ij over j, optionally skipping one input.
inline double prodH(const Net& m, int i, const arma::vec& x, int skip = -1) {
  double p = 1.0;
  for (int j = 0; j < m.n; ++j) {
    if (j == skip || m.reg(i, j) == 0) continue;
    p *= H(m, i, j, x(j));
  }
  return p;
}

arma::vec rhs(const Net& m, const arma::vec& x) {
  arma::vec f(m.n);
  for (int i = 0; i < m.n; ++i)
    f(i) = m.alpha(i) * prodH(m, i, x) - m.gamma(i) * x(i);
  return f;
}

arma::mat jac(const Net& m, const arma::vec& x) {
  arma::mat J(m.n, m.n, arma::fill::zeros);
  for (int i = 0; i < m.n; ++i) {
    for (int k = 0; k < m.n; ++k) {
      int r = m.reg(i, k);
      if (r == 0) continue;
      double d = dact_dx(x(k), m.theta(i, k), m.hill(i, k));
      if (r < 0) d = -d;
      J(i, k) = m.alpha(i) * d * prodH(m, i, x, k);
    }
    J(i, i) -= m.gamma(i);
  }
  return J;
}

arma::mat dfdp(const Net& m, const Loci& loci, const arma::vec& x) {
  arma::mat D(m.n, loci.size(), arma::fill::zeros);
  for (int l = 0; l < loci.size(); ++l) {
    int i = loci.i[l], j = loci.j[l];
    switch (loci.code[l]) {
    case 1: D(i, l) = prodH(m, i, x); break;
    case 2: D(i, l) = -x(i); break;
    case 3:
    case 4: {
      int r = m.reg(i, j);
      if (r == 0) break;
      double d = loci.code[l] == 3 ? dact_dth(x(j), m.theta(i, j), m.hill(i, j))
                                   : dact_dh(x(j), m.theta(i, j), m.hill(i, j));
      if (r < 0) d = -d;
      D(i, l) = m.alpha(i) * d * prodH(m, i, x, j);
      break;
    }
    }
  }
  return D;
}

inline double fnorm(const arma::vec& f, const arma::vec& x) {
  return arma::norm(f, "inf") / (1.0 + arma::norm(x, "inf"));
}

// Dormand-Prince 5(4) integration until the relative residual drops below
// `thresh`. Returns 0 on success, 1 if tmax/maxsteps exceeded, 3 if non-finite.
int integrate(const Net& m, arma::vec& x, double thresh, double tmax,
              int maxsteps, double& dt) {
  static const double a21 = 1.0 / 5, a31 = 3.0 / 40, a32 = 9.0 / 40,
                      a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9,
                      a51 = 19372.0 / 6561, a52 = -25360.0 / 2187,
                      a53 = 64448.0 / 6561, a54 = -212.0 / 729,
                      a61 = 9017.0 / 3168, a62 = -355.0 / 33,
                      a63 = 46732.0 / 5247, a64 = 49.0 / 176,
                      a65 = -5103.0 / 18656, b1 = 35.0 / 384,
                      b3 = 500.0 / 1113, b4 = 125.0 / 192,
                      b5 = -2187.0 / 6784, b6 = 11.0 / 84,
                      e1 = 71.0 / 57600, e3 = -71.0 / 16695,
                      e4 = 71.0 / 1920, e5 = -17253.0 / 339200,
                      e6 = 22.0 / 525, e7 = -1.0 / 40;
  const double atol = 1e-9, rtol = 1e-7;
  double t = 0.0;
  arma::vec k1 = rhs(m, x);
  for (int step = 0; step < maxsteps; ++step) {
    double r = fnorm(k1, x);
    if (!std::isfinite(r)) return 3;
    if (r < thresh) return 0;
    if (t > tmax) return 1;
    arma::vec k2 = rhs(m, x + dt * (a21 * k1));
    arma::vec k3 = rhs(m, x + dt * (a31 * k1 + a32 * k2));
    arma::vec k4 = rhs(m, x + dt * (a41 * k1 + a42 * k2 + a43 * k3));
    arma::vec k5 = rhs(m, x + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
    arma::vec k6 = rhs(m, x + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));
    arma::vec xn = x + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6);
    arma::vec k7 = rhs(m, xn);
    arma::vec err = dt * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7);
    double en = 0.0;
    for (int i = 0; i < m.n; ++i) {
      double sc = atol + rtol * std::max(std::abs(x(i)), std::abs(xn(i)));
      en = std::max(en, std::abs(err(i)) / sc);
    }
    if (!std::isfinite(en)) return 3;
    if (en <= 1.0) {
      t += dt;
      x = arma::clamp(xn, 0.0, arma::datum::inf);
      k1 = rhs(m, x);
    }
    double fac = en > 0 ? 0.9 * std::pow(en, -0.2) : 5.0;
    dt *= std::min(5.0, std::max(0.2, fac));
  }
  return 1;
}

inline bool is_stable(const Net& m, const arma::vec& x) {
  arma::cx_vec ev = arma::eig_gen(jac(m, x));
  return arma::max(arma::real(ev)) < 0;
}

// Newton start: solve f(x) = 0 directly from the starting state x0, for starts
// close to the steady state (e.g. an offspring starting from its parent's
// expression). Accepted only if Newton converges within a few iterations to a
// non-negative, stable steady state within `radius` (relative, max-norm) of
// x0, and the flow at x0 points toward it. The last condition keeps
// development on the attractor that integration from x0 reaches: it rejects
// a jump across a nearby unstable state into another basin. On success x is
// replaced by the steady state and true is returned.
bool newton_start(const Net& m, arma::vec& x, double tol, double radius) {
  arma::vec f0 = rhs(m, x);
  if (!f0.is_finite()) return false;
  arma::vec xn = x, f = f0;
  double lim = radius * (1.0 + arma::norm(x, "inf"));
  bool ok = false;
  for (int it = 0; it < 8; ++it) {
    if (fnorm(f, xn) < tol) { ok = true; break; }
    arma::vec dx;
    if (!arma::solve(dx, jac(m, xn), -f, arma::solve_opts::no_approx)) return false;
    xn += dx;
    if (arma::norm(xn - x, "inf") > lim) return false;
    f = rhs(m, xn);
    if (!f.is_finite()) return false;
  }
  if (!ok || arma::any(xn < -1e-12)) return false;
  arma::vec d = xn - x;
  if (arma::dot(f0, d) < 0) return false;
  arma::vec xc = arma::clamp(xn, 0.0, arma::datum::inf);
  if (!is_stable(m, xc)) return false;
  x = xc;
  return true;
}

// Develop from x (in place) to a steady state. With newton_radius > 0, first
// try a Newton start (see above); otherwise, or if it is not accepted,
// integrate to a small residual and polish with Newton.
// Status: 0 = stable steady state, 1 = not converged, 2 = unstable steady
// state, 3 = non-finite values.
int develop_one(const Net& m, arma::vec& x, double tol, double tmax,
                int maxsteps, bool check_stability, double newton_radius) {
  if (newton_radius > 0 && newton_start(m, x, tol, newton_radius)) return 0;
  double dt = 0.1;
  int s = integrate(m, x, 1e-5, tmax, maxsteps, dt);
  if (s != 0) return s;
  // Newton polish, accepted only if it stays on the same branch.
  arma::vec xn = x, f = rhs(m, x);
  bool ok = false;
  for (int it = 0; it < 30; ++it) {
    arma::vec dx;
    if (!arma::solve(dx, jac(m, xn), -f, arma::solve_opts::no_approx)) break;
    xn += dx;
    f = rhs(m, xn);
    if (!f.is_finite()) break;
    if (fnorm(f, xn) < tol) { ok = true; break; }
  }
  if (ok && arma::all(xn >= -1e-12) &&
      arma::norm(xn - x, "inf") < 1e-3 * (1.0 + arma::norm(x, "inf"))) {
    x = arma::clamp(xn, 0.0, arma::datum::inf);
  } else {
    s = integrate(m, x, tol, tmax, maxsteps, dt);
    if (s != 0) return s;
  }
  if (check_stability && !is_stable(m, x)) return 2;
  return 0;
}

} // namespace

// [[Rcpp::export]]
List hill_develop_cpp(List base, IntegerVector code, IntegerVector li,
                      IntegerVector lj, NumericMatrix G, NumericMatrix X0,
                      double tol, double tmax, int maxsteps,
                      bool check_stability, double newton_radius) {
  Net proto = make_net(base);
  Loci loci{code, li, lj};
  int N = G.nrow(), n = proto.n;
  if (G.ncol() != loci.size()) stop("genotype matrix has wrong number of columns");
  if (X0.nrow() != N || X0.ncol() != n) stop("X0 has wrong dimensions");
  NumericMatrix X(N, n);
  IntegerVector status(N);
  for (int r = 0; r < N; ++r) {
    Net m = proto;
    apply_genotype(m, loci, &G(r, 0), N);
    arma::vec x(n);
    for (int i = 0; i < n; ++i) x(i) = std::max(0.0, X0(r, i));
    status[r] = develop_one(m, x, tol, tmax, maxsteps, check_stability, newton_radius);
    for (int i = 0; i < n; ++i) X(r, i) = x(i);
  }
  return List::create(_["X"] = X, _["status"] = status);
}

// [[Rcpp::export]]
List hill_derivs_cpp(List base, IntegerVector code, IntegerVector li,
                     IntegerVector lj, NumericVector g, NumericVector x) {
  Net m = make_net(base);
  Loci loci{code, li, lj};
  apply_genotype(m, loci, g.begin(), 1);
  arma::vec xv = as<arma::vec>(x);
  return List::create(_["f"] = rhs(m, xv), _["J"] = jac(m, xv),
                      _["dfdp"] = dfdp(m, loci, xv));
}

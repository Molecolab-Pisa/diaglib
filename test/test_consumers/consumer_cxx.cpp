// C++ program using an installed DiagLib through find_package(diaglib) and its C interface.
// Generalized symmetric problem: A = H diag(d) H, B = H diag(s) H -> lambda = d/s.
// Linear response (SMO-GD): A+B = H diag(a) H, A-B = H diag(b) H, S = I, D = 0 -> omega = sqrt(a b).
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>
#include "diaglib.h"

namespace {
constexpr dgl_int N = 800, N_TARG = 3, N_MAX = 6;
std::vector<double> v(N), d(N), s(N), a(N), b(N);
thread_local double scale = 1.0;
int failures = 0;

double dot(const double* x, const double* y) { double r = 0; for (dgl_int i = 0; i < N; ++i) r += x[i] * y[i]; return r; }

void hdh(const std::vector<double>& diag, dgl_int* n, dgl_int* m, double* x, double* y) {
  for (dgl_int j = 0; j < *m; ++j) {
    double* xj = x + j * *n; double* yj = y + j * *n;
    double c = 2 * dot(v.data(), xj);
    for (dgl_int i = 0; i < N; ++i) yj[i] = diag[i] * (xj[i] - c * v[i]);
    c = 2 * dot(v.data(), yj);
    for (dgl_int i = 0; i < N; ++i) yj[i] = scale * (yj[i] - c * v[i]);
  }
}
extern "C" void mv_a(dgl_int* n, dgl_int* m, double* x, double* y) { hdh(d, n, m, x, y); }
extern "C" void mv_b(dgl_int* n, dgl_int* m, double* x, double* y) {
  const double sc = scale; scale = 1.0; hdh(s, n, m, x, y); scale = sc;
}
extern "C" void pc(dgl_int* n, dgl_int* m, double* shift, double* x, double* y) {
  for (dgl_int j = 0; j < *m; ++j)
    for (dgl_int i = 0; i < N; ++i) {
      double den = scale * d[i] / s[i] + *shift;
      if (std::fabs(den) < 1e-3) den = std::copysign(1e-3, den);
      y[i + j * *n] = x[i + j * *n] / den;
    }
}
extern "C" void apb(dgl_int* n, dgl_int* m, double* x, double* y) { hdh(a, n, m, x, y); }
extern "C" void amb(dgl_int* n, dgl_int* m, double* x, double* y) { hdh(b, n, m, x, y); }
extern "C" void ident(dgl_int* n, dgl_int* m, double* x, double* y) { std::memcpy(y, x, sizeof(double) * *n * *m); }
extern "C" void lrprec(dgl_int* n, dgl_int* m, double* fac, double* xp, double* xm, double* yp, double* ym) {
  const double f = *fac;
  for (dgl_int j = 0; j < *m; ++j)
    for (dgl_int i = 0; i < N; ++i) {
      double den = f * f * a[i] * b[i] - 1.0;
      if (std::fabs(den) < 1e-3) den = 1e-3;
      yp[i + j * *n] = (f * b[i] * xp[i + j * *n] + xm[i + j * *n]) / den;
      ym[i + j * *n] = (f * a[i] * xm[i + j * *n] + xp[i + j * *n]) / den;
    }
}

std::vector<double> lowest(std::vector<double> x, double factor) {
  for (auto& e : x) e *= factor;
  std::sort(x.begin(), x.end());
  return {x.begin(), x.begin() + N_TARG};
}
// relative residual of a paired linear response equation t = omega x
double pair_err(const std::vector<double>& t, const std::vector<double>& x, double omega) {
  double r = 0;
  for (dgl_int i = 0; i < N; ++i) r += std::pow(t[i] - omega * x[i], 2);
  return std::sqrt(r) / (std::sqrt(dot(t.data(), t.data())) + omega * std::sqrt(dot(x.data(), x.data())));
}
double max_err(const std::vector<double>& eig, const std::vector<double>& ref) {
  double e = 0; for (dgl_int k = 0; k < N_TARG; ++k) e = std::max(e, std::fabs(eig[k] - ref[k])); return e;
}
void check(const char* label, bool cond) {
  if (!cond) ++failures;
  std::printf("%-45s %s\n", label, cond ? "PASSED" : "FAILED");
}
}  // namespace

int main() {
  double nv = 0;
  for (dgl_int i = 0; i < N; ++i) {
    v[i] = std::sin(i + 1.0); nv += v[i] * v[i];
    d[i] = 1.0 + i + 0.5 * std::pow(std::sin(i + 1.0), 2);
    s[i] = 1.0 + 0.2 * std::pow(std::cos(i + 1.0), 2);
    a[i] = 2.0 + i; b[i] = 1.0 + 0.5 * (i + 1);
  }
  for (auto& x : v) x /= std::sqrt(nv);
  std::printf("header integer kind %d, library %d\n", DGL_INT_KIND, dgl_integer_kind());

  std::vector<double> eig(N_MAX), evec(N * N_MAX, 0.0), ref_gen(N);
  bool ok = false; dgl_int info = -1;
  for (dgl_int i = 0; i < N; ++i) ref_gen[i] = d[i] / s[i];
  dgl_davidson_driver(N, N_TARG, N_MAX, mv_a, pc, mv_b, eig.data(), evec.data(), &ok, &info, false, 1e-9, 100, 20, 0.0, true, 50, "MB");
  std::printf("  eig err %.2e\n", max_err(eig, lowest(ref_gen, 1.0)));
  check("generalized Davidson", ok && info == DGL_SUCCESS && max_err(eig, lowest(ref_gen, 1.0)) < 1e-8);

  std::vector<double> evec2(2 * N * N_MAX, 0.0), ref_lr(N);
  for (dgl_int i = 0; i < N; ++i) ref_lr[i] = std::sqrt(a[i] * b[i]);
  dgl_smogd_driver(2 * N, N_TARG, N_MAX, apb, amb, ident, ident, lrprec, eig.data(), evec2.data(), &ok, &info, false, 1e-9, 200, 20, 50, "MB");
  // the whole problem, [A B; B A] (y;z) = omega [1 0; 0 -1] (y;z), is the pair
  // (A+B) x+ = omega x-, (A-B) x- = omega x+ with x+ = y + z and x- = y - z: the sum and the
  // difference of its two block rows. That is the form the callbacks give (A and B are only
  // available through A+B and A-B) and the one the driver converges. Eliminating x+ instead
  // gives (A+B)(A-B) x- = omega^2 x-, whose residual is these ones amplified by |A+B|/omega^2,
  // about 270 here: too strict to be compared with the requested tolerance.
  std::vector<double> xp(N), xm(N), t(N); dgl_int n = N, one = 1;
  double res = 0;
  for (dgl_int k = 0; k < N_TARG; ++k) {
    for (dgl_int i = 0; i < N; ++i) {
      xp[i] = evec2[i + k * 2 * N] + evec2[N + i + k * 2 * N];
      xm[i] = evec2[i + k * 2 * N] - evec2[N + i + k * 2 * N];
    }
    apb(&n, &one, xp.data(), t.data());
    res = std::max(res, pair_err(t, xm, eig[k]));
    amb(&n, &one, xm.data(), t.data());
    res = std::max(res, pair_err(t, xp, eig[k]));
  }
  std::printf("  eig err %.2e paired residual %.2e\n", max_err(eig, lowest(ref_lr, 1.0)), res);
  check("SMO-GD", ok && info == DGL_SUCCESS && max_err(eig, lowest(ref_lr, 1.0)) < 1e-8 && res < 1e-6);

  dgl_smogd_driver(2 * N - 1, N_TARG, N_MAX, apb, amb, ident, ident, lrprec, eig.data(), evec2.data(), &ok, &info, false, 1e-9, 200, 20, 50, "MB");
  check("error: odd size for SMO-GD", !ok && info == DGL_ERR_INPUT);

  // std::thread: generalized Davidson with different scales at the same time
  std::vector<std::thread> threads;
  std::atomic<int> passed{0};
  std::vector<double> errs(8);
  for (int t = 0; t < 8; ++t)
    threads.emplace_back([t, &passed, &errs, &ref_gen] {
      scale = 1.0 + t;
      std::vector<double> e(N_MAX), ev(N * N_MAX, 0.0);
      bool lok = false; dgl_int linfo = -1;
      dgl_davidson_driver(N, N_TARG, N_MAX, mv_a, pc, mv_b, e.data(), ev.data(), &lok, &linfo, false, 1e-9, 400, 20, 0.0, true, 50, "MB");
      errs[t] = max_err(e, lowest(ref_gen, scale));
      if (lok && linfo == DGL_SUCCESS && errs[t] < 1e-8) ++passed;
      else std::printf("  FAILED thread %d: ok=%d info=%ld eigenvalue error=%.2e\n", t, lok, (long)linfo, errs[t]);
    });
  for (auto& th : threads) th.join();
  std::printf("  max eigenvalue error %.2e\n", *std::max_element(errs.begin(), errs.end()));
  check("8 concurrent std::threads", passed == 8);

  std::printf("C++ consumer: %d failures\n", failures);
  return failures != 0;
}

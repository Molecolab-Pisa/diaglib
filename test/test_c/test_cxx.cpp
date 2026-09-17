//
// check that diaglib.h can be included and the library linked from C++:
// solve a small symmetric problem with Davidson.
//
#include <cmath>
#include <cstdio>
#include <vector>
#include "diaglib.h"

static void matvec(dgl_int* n, dgl_int* m, double* x, double* ax) {
  const dgl_int N = *n, M = *m;
  for (dgl_int j = 0; j < M; ++j)
    for (dgl_int i = 0; i < N; ++i) {
      double s = 0.0;
      for (dgl_int k = 0; k < N; ++k)
        s += ((i == k) ? (i + 2.0) : 1.0 / (i + k + 2.0)) * x[k + j * N];
      ax[i + j * N] = s;
    }
}

static void precnd(dgl_int* n, dgl_int* m, double* shift, double* r, double* z) {
  const dgl_int N = *n, M = *m;
  for (dgl_int j = 0; j < M; ++j)
    for (dgl_int i = 0; i < N; ++i) z[i + j * N] = r[i + j * N] / (i + 2.0 + *shift);
}

int main() {
  const dgl_int n = 200, n_targ = 2, n_max = 4;
  std::vector<double> eig(n_max), evec(n * n_max, 0.0);
  for (dgl_int i = 0; i < n_max; ++i) evec[i + i * n] = 1.0;
  bool ok = false;
  dgl_int info = 0;
  dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, nullptr, eig.data(), evec.data(), &ok, &info,
                      false, 1e-8, 100, 25, 0.0, 10, "MB");
  // reference value computed with LAPACK
  const bool pass = ok && info == DGL_SUCCESS && std::fabs(eig[0] - 1.86940073) < 1e-7;
  std::printf("C++ test: ok=%d info=%ld eig=%.8f -> %s\n", ok, static_cast<long>(info), eig[0], pass ? "PASSED" : "FAILED");
  return pass ? 0 : 1;
}

#include <stdio.h>
#include <stdbool.h>
#include "diaglib.h"

void matvec_c(int* n, int* m, double* x, double* ax) {
  int N = *n;
  int M = *m;

  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int k = 0; k < N; ++k) {
        double A_ik = (i == k) ? (i + 2.0) : (1.0 / (i + k + 2.0)); // i,k zero-based
        sum += A_ik * x[k + j * N];
      }
      ax[i + j * N] = sum;
    }
  }
}

void precnd_c(int* n, int* m, double* shift, double* r, double* z) {
  int N = *n;
  int M = *m;

  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double diag = (i + 1.0) + (*shift); // A_ii + shift
      z[i + j * N] = r[i + j * N] / diag;
    }
  }
}

void bvec_c(int* n, int* m, double* x, double* bx) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < N; ++i)
      bx[i + j * N] = x[i + j * N];  // identità
}

int main() {
  const int n = 1000, n_targ = 20, n_max = 25, max_iter = 100, max_dav = 20;
  const double tol = 1e-8, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max];
  bool ok;

  printf("\nCalling DAVIDSON driver...\n");
  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;

  davidson_driver_c(true, n, n_targ, n_max, max_iter, max_dav, tol, shift,
                    matvec_c, precnd_c, eig, evec, &ok);

  if (ok) {
    printf("Davidson converged.\nEigenvalues:\n");
    for (int i = 0; i < n_targ; ++i)
      printf("  %.12f\n", eig[i]);
  } else {
    printf("Davidson failed to converge.\n");
  }

  printf("\nCalling LOBPCG driver...\n");

  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;

  ok = false;
  
  lobpcg_driver_c(
    true,         // verbose
    false,        // gen_eig = false → problema standard
    n, n_targ, n_max, max_iter, tol, shift,
    matvec_c, precnd_c, bvec_c,
    eig, evec, &ok
  );
  
  if (ok) {
    printf("LOBPCG converged.\nEigenvalues:\n");
    for (int i = 0; i < n_targ; ++i)
      printf("  %.12f\n", eig[i]);
  } else {
    printf("LOBPCG failed to converge.\n");
  }

  return 0;
}
